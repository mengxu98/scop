import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile


def git(*args):
    return subprocess.check_output(["git", *args], text=True).strip()


def fetch(ref):
    return subprocess.run(
        ["git", "fetch", "--no-tags", "--depth=1", "origin", ref],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0


def pages(root):
    topics = sorted(p.stem for p in (root / "man").glob("*.Rd"))
    articles = sorted(
        p.relative_to(root / "vignettes").with_suffix("").as_posix()
        for p in (root / "vignettes").rglob("*.Rmd")
        if not any(part.startswith("_") for part in p.relative_to(root / "vignettes").parts)
    )
    return topics, articles


def plan_changes(root, changes, full=False):
    all_topics, all_articles = pages(root)
    topics, articles, deleted = set(), set(), set()
    home = news = reference_index = articles_index = False
    reasons = []
    for status, path in changes:
        if path.startswith(("R/", "src/", "data/", "inst/")) or path in (
            "DESCRIPTION", "NAMESPACE", "renv.lock", "renv/settings.json", "_pkgdown.yml"
        ) or path.startswith(("pkgdown/", ".github/pkgdown-")) or path == ".github/workflows/pkgdown.yaml":
            full = True
            reasons.append(path)
        elif path.startswith("man/") and path.endswith(".Rd"):
            topic = PurePosixPath(path).stem
            reference_index = True
            if status == "D":
                deleted.add(f"reference/{topic}.html")
            else:
                topics.add(topic)
        elif path.startswith("vignettes/"):
            articles_index = True
            if path.endswith(".Rmd"):
                article = str(PurePosixPath(path).relative_to("vignettes").with_suffix(""))
                if any(part.startswith("_") for part in PurePosixPath(article).parts):
                    articles.update(all_articles)
                elif status == "D":
                    deleted.add(f"articles/{article}.html")
                else:
                    articles.add(article)
            else:
                articles.update(all_articles)
        elif path.startswith("README.") or path.startswith("man/figures/"):
            home = True
        elif path == "NEWS.md":
            news = True
    if full:
        topics, articles = set(all_topics), set(all_articles)
        home = news = reference_index = articles_index = True
    execute = bool(topics or articles)
    return {
        "full": full,
        "topics": sorted(topics),
        "articles": sorted(articles),
        "home": home,
        "news": news,
        "reference_index": reference_index,
        "articles_index": articles_index,
        "deleted": sorted(deleted),
        "execute": execute,
        "build": bool(execute or home or news or reference_index or articles_index or deleted),
        "reasons": reasons,
    }


def restore(plan):
    destination = Path("docs")
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir()
    if plan["full"]:
        return
    if not fetch("gh-pages") or git("rev-parse", "FETCH_HEAD") != plan["published"]:
        raise RuntimeError("Published website changed after planning; rerun against its current revision")
    with tempfile.TemporaryFile() as archive:
        subprocess.run(["git", "archive", "FETCH_HEAD"], stdout=archive, check=True)
        archive.seek(0)
        with tarfile.open(fileobj=archive) as tree:
            tree.extractall(destination, filter="data")
    for name in plan["deleted"]:
        (destination / name).unlink(missing_ok=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "restore":
        restore(json.loads(Path("pkgdown-plan.json").read_text()))
        return
    if len(sys.argv) > 1 and sys.argv[1] == "full":
        plan = json.loads(Path("pkgdown-plan.json").read_text())
        deleted, reasons = plan["deleted"], plan["reasons"]
        plan.update(plan_changes(Path("."), [], full=True))
        plan.update(deleted=deleted, reasons=reasons + ["R or installed dependency versions changed"])
        Path("pkgdown-plan.json").write_text(json.dumps(plan, indent=2) + "\n")
        with open(os.environ["GITHUB_ENV"], "a") as stream:
            for key in ("topics", "articles"):
                stream.write(f"PKGDOWN_{key.upper()}={','.join(plan[key])}\n")
            stream.write("PKGDOWN_EXECUTE=true\n")
        return
    published = None
    manifest = {}
    if fetch("gh-pages"):
        published = git("rev-parse", "FETCH_HEAD")
        try:
            manifest = json.loads(git("show", "FETCH_HEAD:.pkgdown-build.json"))
        except (subprocess.CalledProcessError, ValueError):
            pass
    event = os.environ.get("PKGDOWN_EVENT", "workflow_dispatch")
    base = manifest.get("source")
    full = os.environ.get("PKGDOWN_FORCE") == "true" or event in ("schedule", "release")
    changes = []
    if not base or not fetch(base):
        full = True
    else:
        lines = git("diff", "--no-renames", "--name-status", base, "HEAD").splitlines()
        changes = [line.split("\t", 1) for line in lines]
    if not manifest:
        full = True
    plan = plan_changes(Path("."), changes, full)
    plan.update(source=git("rev-parse", "HEAD"), published=published, base=base)
    Path("pkgdown-plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    print(json.dumps(plan, indent=2))
    output = os.environ.get("GITHUB_OUTPUT")
    if output:
        with open(output, "a") as stream:
            for key in ("build", "execute"):
                stream.write(f"{key}={str(plan[key]).lower()}\n")
            for key in ("topics", "articles"):
                stream.write(f"{key}={','.join(plan[key])}\n")
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as stream:
            stream.write(f"Website build: {'full' if plan['full'] else 'incremental'}\n\n")
            stream.write(f"Reference topics: {len(plan['topics'])}; articles: {len(plan['articles'])}\n")


if __name__ == "__main__":
    main()
