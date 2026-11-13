import importlib.util
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


def load_runner(path):
    spec = importlib.util.spec_from_file_location("scop_official_scenicplus", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def wait_for(predicate, label, seconds=30):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.05)
    raise AssertionError("Timed out waiting for " + label)


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    state = subprocess.run(["ps", "-p", str(pid), "-o", "stat="], capture_output=True, text=True, check=False).stdout.strip()
    return bool(state) and not state.startswith("Z")


def child_mode(runner_path, mode, directory):
    module = load_runner(runner_path)
    directory = Path(directory)
    if mode == "try-lock":
        try:
            with module.execution_lock(directory):
                module.atomic_json(directory / "acquired.json", dict(pid=os.getpid()))
        except RuntimeError:
            raise SystemExit(17)
        return
    if mode == "hold":
        with module.execution_lock(directory):
            module.atomic_json(directory / "ready.json", dict(pid=os.getpid()))
            wait_for(lambda: (directory / "release").is_file(), "release", 60)
        return
    if mode == "orphan":
        with module.execution_lock(directory) as lock_fd:
            code = "import json,os,time;from pathlib import Path;Path(" + repr(str(directory / "orphan-child.json")) + ").write_text(json.dumps({'pid':os.getpid()}));time.sleep(60)"
            process = subprocess.Popen([sys.executable, "-c", code], pass_fds=(lock_fd,), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            wait_for(lambda: (directory / "orphan-child.json").exists(), "orphan child")
            module.atomic_json(directory / "orphan-parent.json", dict(pid=os.getpid(), child=process.pid))
            while True:
                time.sleep(1)
    if mode == "workflow":
        with module.execution_lock(directory) as lock_fd:
            code = "import json,os,subprocess,sys,time;from pathlib import Path;p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)']);Path(" + repr(str(directory / "workflow-child.json")) + ").write_text(json.dumps({'pid':os.getpid(),'child':p.pid}));time.sleep(60);p.wait()"
            before = {sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)}
            try:
                module.run_workflow([sys.executable, "-c", code], directory, os.environ.copy(), lock_fd)
            except module.WorkflowCancelled as error:
                assert all(signal.getsignal(sig) == handler for sig, handler in before.items())
                module.atomic_json(directory / "cancelled.json", dict(signal=error.signum, status="cancelled"))
                raise SystemExit(128 + error.signum)
        return
    if mode == "workflow-early-exit":
        with module.execution_lock(directory) as lock_fd:
            child_code = "import json,os,signal,time;from pathlib import Path;signal.signal(signal.SIGTERM,signal.SIG_IGN);Path(" + repr(str(directory / "grandchild.json")) + ").write_text(json.dumps({'pid':os.getpid()}));time.sleep(60)"
            code = "import json,os,subprocess,sys,time;from pathlib import Path;p=subprocess.Popen([sys.executable,'-c'," + repr(child_code) + "]);f=Path(" + repr(str(directory / "grandchild.json")) + ");deadline=time.monotonic()+10\nwhile not f.exists() and time.monotonic()<deadline:time.sleep(.02)\nassert f.exists()"
            status = module.run_workflow([sys.executable, "-c", code], directory, os.environ.copy(), lock_fd)
            assert status == 0
        return
    if mode == "workflow-clean-exit":
        with module.execution_lock(directory) as lock_fd:
            status = module.run_workflow([sys.executable, "-c", "raise SystemExit(0)"], directory, os.environ.copy(), lock_fd)
            assert status == 0
        return
    raise ValueError("Unknown test mode")


def probe(runner_path, directory, mode, wait=True):
    command = [sys.executable, str(Path(__file__).resolve()), runner_path, mode, str(directory)]
    if wait:
        return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    return subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main(runner_path, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    module = load_runner(runner_path)
    if os.name != "posix":
        try:
            with module.execution_lock(output):
                raise AssertionError("Unsupported platform acquired a POSIX lock")
        except RuntimeError as error:
            assert "Linux or macOS" in str(error)
        print("UNSUPPORTED PLATFORM MESSAGE PASSED")
        return
    records = []
    directory = output / "concurrent"
    directory.mkdir()
    sentinel = directory / "request_identity.json"
    sentinel.write_text('{"existing":"unchanged"}')
    original = sentinel.read_bytes()
    worker = probe(runner_path, directory, "hold", wait=False)
    try:
        wait_for(lambda: (directory / "ready.json").exists(), "first lock owner")
        inode = (directory / ".scop-execution.lock").stat().st_ino
        second = probe(runner_path, directory, "try-lock")
        assert second.returncode == 17 and not (directory / "acquired.json").exists()
        assert sentinel.read_bytes() == original
        (directory / "release").touch()
        assert worker.wait(timeout=30) == 0
        third = probe(runner_path, directory, "try-lock")
        assert third.returncode == 0, third.stderr
        assert (directory / ".scop-execution.lock").stat().st_ino == inode
        records.append(dict(case="different_process_rejected_before_write_then_same_inode_reacquired", passed=True))
    finally:
        if worker.poll() is None:
            worker.kill()
            worker.wait()
    for signum in (signal.SIGINT, signal.SIGTERM):
        directory = output / signal.Signals(signum).name
        directory.mkdir()
        worker = probe(runner_path, directory, "workflow", wait=False)
        pids = []
        try:
            wait_for(lambda: (directory / "workflow-child.json").exists(), "real workflow subprocess")
            info = json.loads((directory / "workflow-child.json").read_text())
            pids = [info["pid"], info["child"]]
            os.kill(worker.pid, signum)
            assert worker.wait(timeout=30) == 128 + signum
            wait_for(lambda: not any(alive(pid) for pid in pids), "workflow process group to exit")
            receipt = json.loads((directory / "cancelled.json").read_text())
            assert receipt == dict(signal=signum, status="cancelled")
            second = probe(runner_path, directory, "try-lock")
            assert second.returncode == 0, second.stderr
            records.append(dict(case=signal.Signals(signum).name, passed=True, child_and_grandchild_exited=True, handlers_restored=True))
        finally:
            if worker.poll() is None:
                worker.kill()
                worker.wait()
            for pid in pids:
                if alive(pid):
                    os.kill(pid, signal.SIGKILL)
    directory = output / "early-exit"
    directory.mkdir()
    worker = probe(runner_path, directory, "workflow-early-exit", wait=False)
    child = None
    try:
        wait_for(lambda: (directory / "grandchild.json").exists(), "TERM-ignoring grandchild")
        child = json.loads((directory / "grandchild.json").read_text())["pid"]
        assert worker.wait(timeout=30) == 0
        wait_for(lambda: not alive(child), "grandchild after early leader exit", seconds=5)
        second = probe(runner_path, directory, "try-lock")
        assert second.returncode == 0, second.stderr
        records.append(dict(case="early_leader_exit_cleans_TERM_ignoring_grandchild", passed=True))
    finally:
        if worker.poll() is None:
            worker.kill()
            worker.wait()
        if child is not None and alive(child):
            os.kill(child, signal.SIGKILL)
    directory = output / "clean-exit"
    directory.mkdir()
    clean = probe(runner_path, directory, "workflow-clean-exit")
    assert clean.returncode == 0, clean.stderr
    records.append(dict(case="clean_workflow_exit_preserves_returncode", passed=True))
    directory = output / "orphan"
    directory.mkdir()
    worker = probe(runner_path, directory, "orphan", wait=False)
    child = None
    try:
        wait_for(lambda: (directory / "orphan-parent.json").exists(), "orphan parent readiness")
        info = json.loads((directory / "orphan-parent.json").read_text())
        child = info["child"]
        os.kill(worker.pid, signal.SIGKILL)
        assert worker.wait(timeout=30) == -signal.SIGKILL
        assert alive(child)
        second = probe(runner_path, directory, "try-lock")
        assert second.returncode == 17
        os.kill(child, signal.SIGTERM)
        wait_for(lambda: not alive(child), "orphan child exit")
        third = probe(runner_path, directory, "try-lock")
        assert third.returncode == 0, third.stderr
        records.append(dict(case="SIGKILL_parent_retains_child_lock_until_child_exit", passed=True))
    finally:
        if child is not None and alive(child):
            os.kill(child, signal.SIGKILL)
    path = output / "atomic.json"
    module.atomic_json(path, dict(complete=True, values=list(range(100))))
    original = path.read_bytes()
    try:
        module.atomic_json(path, dict(invalid=object()))
        raise AssertionError("An unserializable value was accepted")
    except TypeError:
        pass
    assert path.read_bytes() == original and not list(output.glob(".atomic.json.*.tmp"))
    module.atomic_json(path, dict(replaced=True))
    assert json.loads(path.read_text()) == dict(replaced=True)
    records.append(dict(case="atomic_json_preserves_previous_value_on_failure_and_cleans_stage", passed=True))
    directory = output / "foreign-directory"
    directory.mkdir()
    foreign = directory / "native-output.txt"
    foreign.write_text("foreign native result")
    original = foreign.read_bytes()
    try:
        with module.execution_lock(directory) as lock_fd:
            module.run_locked({}, {}, directory, directory, lock_fd)
        raise AssertionError("An unverified nonempty workspace was accepted")
    except ValueError as error:
        assert "empty or contain a verified scop request identity" in str(error)
    assert foreign.read_bytes() == original and not (directory / "request_identity.json").exists()
    records.append(dict(case="foreign_nonempty_workspace_rejected_without_mutation", passed=True))
    directory = output / "unverified-native-locks"
    native_locks = directory / ".snakemake" / "locks"
    native_locks.mkdir(parents=True)
    (native_locks / "0.input.lock").write_text("input\n")
    (native_locks / "0.output.lock").write_text("output\n")
    snapshot = module.snakemake_lock_snapshot(directory)
    try:
        module.verify_resume_locks(directory, "identity", True)
        raise AssertionError("Native locks without ownership were accepted")
    except ValueError as error:
        assert "official manual unlock" in str(error)
    assert module.snakemake_lock_snapshot(directory) == snapshot
    records.append(dict(case="unverified_native_locks_rejected_without_mutation", passed=True))
    module.atomic_json(directory / "snakemake_lock_owner.json", dict(identity_sha256="identity", pgid=os.getpgrp(), native_locks=snapshot))
    try:
        module.verify_resume_locks(directory, "identity", True)
        raise AssertionError("A live process group was accepted as stale")
    except ValueError as error:
        assert "official manual unlock" in str(error)
    assert module.snakemake_lock_snapshot(directory) == snapshot
    records.append(dict(case="live_recorded_PGID_conservatively_rejected", passed=True))
    (native_locks / "0.output.lock").write_text("changed output\n")
    changed = module.snakemake_lock_snapshot(directory)
    try:
        module.verify_resume_locks(directory, "identity", True)
        raise AssertionError("Changed native locks were accepted as owned")
    except ValueError as error:
        assert "official manual unlock" in str(error)
    assert module.snakemake_lock_snapshot(directory) == changed
    records.append(dict(case="changed_native_lock_content_and_identity_rejected", passed=True))
    directory = output / "owned-lock-recovery"
    native_locks = directory / ".snakemake" / "locks"
    native_locks.mkdir(parents=True)
    for kind in ("input", "output"):
        (native_locks / ("0." + kind + ".lock")).write_text(kind + "\n")
    terminated = subprocess.Popen([sys.executable, "-c", "pass"], start_new_session=True)
    terminated.wait(timeout=10)
    snapshot = module.snakemake_lock_snapshot(directory)
    module.atomic_json(directory / "snakemake_lock_owner.json", dict(identity_sha256="identity", pgid=terminated.pid, native_locks=snapshot))
    original_verify = module.verify_resume_locks
    def add_foreign_after_check(*args):
        verified = original_verify(*args)
        for kind in ("input", "output"):
            (native_locks / ("1." + kind + ".lock")).write_text("foreign " + kind + "\n")
        return verified
    module.verify_resume_locks = add_foreign_after_check
    try:
        module.clear_owned_snakemake_locks(directory, "identity", snapshot)
    finally:
        module.verify_resume_locks = original_verify
    remaining = module.snakemake_lock_snapshot(directory)
    assert set(remaining) == {"1.input.lock", "1.output.lock"}
    assert (native_locks / "1.input.lock").read_text() == "foreign input\n"
    records.append(dict(case="foreign_lock_added_after_final_check_is_preserved", passed=True))
    directory = output / "replaced-owned-lock"
    native_locks = directory / ".snakemake" / "locks"
    native_locks.mkdir(parents=True)
    for kind in ("input", "output"):
        (native_locks / ("0." + kind + ".lock")).write_text(kind + "\n")
    snapshot = module.snakemake_lock_snapshot(directory)
    module.atomic_json(directory / "snakemake_lock_owner.json", dict(identity_sha256="identity", pgid=terminated.pid, native_locks=snapshot))
    replaced = native_locks / "0.output.lock"
    staged = directory / "replacement"
    staged.write_bytes(replaced.read_bytes())
    os.utime(staged, ns=(replaced.stat().st_atime_ns, replaced.stat().st_mtime_ns))
    os.replace(staged, replaced)
    changed = module.snakemake_lock_snapshot(directory)
    assert changed[replaced.name]["st_ino"] != snapshot[replaced.name]["st_ino"]
    try:
        module.clear_owned_snakemake_locks(directory, "identity", snapshot)
        raise AssertionError("A replacement native inode was removed")
    except ValueError:
        pass
    assert module.snakemake_lock_snapshot(directory) == changed
    records.append(dict(case="same_name_size_mtime_replacement_inode_rejected_before_removal", passed=True))
    directory = output / "exception"
    directory.mkdir()
    try:
        with module.execution_lock(directory):
            raise ValueError("expected failure")
    except ValueError:
        pass
    with module.execution_lock(directory):
        records.append(dict(case="exception_releases_directory_lock", passed=True))
    result = dict(status="passed", runner=str(Path(runner_path).resolve()), cases=records)
    (output / "receipt.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    if len(sys.argv) == 4:
        child_mode(*sys.argv[1:])
    else:
        main(*sys.argv[1:])
