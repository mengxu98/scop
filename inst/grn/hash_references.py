from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path
import sys


BLOCK_SIZE = 64 * 1024**2
READ_SIZE = 2 * 1024**2


def identity(stat):
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def whole_file_digest(path):
    digest = sha256()
    with path.open("rb", buffering=0) as handle:
        initial = identity(os.fstat(handle.fileno()))
        for data in iter(lambda: handle.read(READ_SIZE), b""):
            digest.update(data)
    if identity(path.stat()) != initial:
        raise OSError(f"Resource changed during verification: {path}")
    return digest.hexdigest()


def file_digest(path, pool):
    state = path.stat()
    initial = identity(state)

    def digest_block(offset):
        count = min(BLOCK_SIZE, state.st_size - offset)
        digest = sha256()
        with path.open("rb", buffering=0) as handle:
            if identity(os.fstat(handle.fileno())) != initial:
                raise OSError(f"Resource changed during verification: {path}")
            handle.seek(offset)
            while count:
                data = handle.read(min(count, READ_SIZE))
                if not data:
                    raise OSError(f"Resource ended during verification: {path}")
                digest.update(data)
                count -= len(data)
        return digest.digest()

    digest = sha256(b"SCOP file identity v1\0" + state.st_size.to_bytes(8, "big") + BLOCK_SIZE.to_bytes(8, "big"))
    for block in pool.map(digest_block, range(0, state.st_size, BLOCK_SIZE)):
        digest.update(block)
    if identity(path.stat()) != initial:
        raise OSError(f"Resource changed during verification: {path}")
    return digest.hexdigest()


def main(algorithm, paths):
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=min(8, cpus)) as pool:
        if algorithm == "sha256":
            result = list(pool.map(whole_file_digest, map(Path, paths)))
        elif algorithm == "sha256-tree-64m-v1":
            result = [file_digest(Path(path), pool) for path in paths]
        else:
            raise ValueError(f"Unsupported resource fingerprint algorithm: {algorithm}")
    print(json.dumps(result))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
