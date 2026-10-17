import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location("planner", Path(__file__).with_name("pkgdown-plan.py"))
planner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(planner)


class WebsitePlanTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        for name in ("man/One.Rd", "man/Two.Rd", "vignettes/first.Rmd", "vignettes/second.Rmd"):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()

    def plan(self, *changes):
        return planner.plan_changes(self.root, changes)

    def test_reference_change_runs_only_affected_examples(self):
        plan = self.plan(("M", "man/One.Rd"))
        self.assertEqual(plan["topics"], ["One"])
        self.assertEqual(plan["articles"], [])
        self.assertTrue(plan["execute"])
        self.assertFalse(plan["full"])

    def test_home_and_reference_changes_are_both_retained(self):
        plan = self.plan(("M", "README.md"), ("M", "man/Two.Rd"))
        self.assertTrue(plan["home"])
        self.assertEqual(plan["topics"], ["Two"])

    def test_news_does_not_require_package_compilation(self):
        plan = self.plan(("M", "NEWS.md"))
        self.assertTrue(plan["build"])
        self.assertTrue(plan["news"])
        self.assertFalse(plan["execute"])

    def test_source_data_and_environment_changes_invalidate_results(self):
        for name in ("R/helper.R", "src/helper.cpp", "data/example.rda", "renv.lock", "DESCRIPTION", "_pkgdown.yml"):
            with self.subTest(name=name):
                plan = self.plan(("M", name))
                self.assertTrue(plan["full"])
                self.assertEqual(plan["topics"], ["One", "Two"])
                self.assertEqual(plan["articles"], ["first", "second"])

    def test_article_change_runs_only_that_article(self):
        plan = self.plan(("M", "vignettes/first.Rmd"))
        self.assertEqual(plan["articles"], ["first"])
        self.assertEqual(plan["topics"], [])

    def test_shared_article_resources_invalidate_articles(self):
        for name in ("vignettes/images/plot.png", "vignettes/_shared.Rmd", "vignettes/data/input.csv"):
            with self.subTest(name=name):
                plan = self.plan(("M", name))
                self.assertEqual(plan["articles"], ["first", "second"])
                self.assertFalse(plan["full"])

    def test_deletions_remove_published_pages(self):
        plan = self.plan(("D", "man/Removed.Rd"), ("D", "vignettes/removed.Rmd"))
        self.assertEqual(plan["deleted"], ["articles/removed.html", "reference/Removed.html"])
        self.assertFalse(plan["execute"])
        self.assertTrue(plan["build"])

    def test_unrelated_changes_skip_r_setup(self):
        self.assertFalse(self.plan(("M", "tests/testthat/test-other.R"))["build"])

    def test_push_uses_last_successful_publication(self):
        original = Path.cwd()
        os.chdir(self.root)
        self.addCleanup(os.chdir, original)

        def git(*args):
            if args == ("rev-parse", "FETCH_HEAD"):
                return "published-tree"
            if args[0] == "show":
                return json.dumps({"source": "last-successful-build"})
            if args[0] == "diff":
                self.assertEqual(args[3], "last-successful-build")
                return "M\tman/One.Rd\nM\tman/Two.Rd"
            return "current-source"

        with patch.object(planner, "fetch", return_value=True), patch.object(planner, "git", side_effect=git), patch.dict(os.environ, {"PKGDOWN_EVENT": "push", "PKGDOWN_FORCE": "false", "GITHUB_OUTPUT": "", "GITHUB_STEP_SUMMARY": ""}):
            planner.main()
        plan = json.loads(Path("pkgdown-plan.json").read_text())
        self.assertEqual(plan["topics"], ["One", "Two"])
        self.assertFalse(plan["full"])


if __name__ == "__main__":
    unittest.main()
