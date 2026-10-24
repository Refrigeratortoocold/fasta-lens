import csv
import gzip
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CliTests(unittest.TestCase):
    def run_cli(self, *args, input=None):
        return subprocess.run([sys.executable, "-m", "fasta_lens", *map(str, args)],
                              input=input, text=True, encoding="utf-8", capture_output=True, cwd=ROOT)

    def test_stdin_json(self):
        result = self.run_cli("-", "--format", "json", input=">x\nACGT\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["summary"]["gc_pct"], 50)
        self.assertEqual(result.stderr, "")

    def test_markdown(self):
        result = self.run_cli("examples/demo.fasta")
        self.assertEqual(result.returncode, 0)
        self.assertIn("Records: 3 | Bases: 30 | N50: 10 bp", result.stdout)
        self.assertIn("above_max_ambiguity", result.stdout)

    def test_tsv(self):
        result = self.run_cli("-", "--format", "tsv", input=">x\nNN\n")
        rows = list(csv.DictReader(io.StringIO(result.stdout), delimiter="\t"))
        self.assertEqual(rows[0]["gc_pct"], "NA")
        self.assertEqual(rows[0]["ambiguity_pct"], "100.00")

    def test_escaped_markdown(self):
        result = self.run_cli("-", input="><img|x>\nAC\n")
        self.assertIn("&lt;img&#124;x&gt;", result.stdout)
        self.assertNotIn("<img", result.stdout)

    def test_fail_on_warning(self):
        result = self.run_cli("-", "--fail-on-warning", "--format", "json", input=">x\nNN\n")
        self.assertEqual(result.returncode, 1)
        self.assertEqual(json.loads(result.stdout)["summary"]["warning_record_count"], 1)

    def test_accepted_threshold(self):
        result = self.run_cli("-", "--fail-on-warning", "--max-ambiguity-pct", "100", input=">x\nNN\n")
        self.assertEqual(result.returncode, 0)

    def test_invalid_input_has_no_partial_report(self):
        result = self.run_cli("-", input=">good\nACGT\n>bad\nU\n")
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, "")
        self.assertIn("line 4", result.stderr)

    def test_invalid_options(self):
        for args in [("--min-length", "0"), ("--max-ambiguity-pct", "nan"), ("--max-ambiguity-pct", "101")]:
            with self.subTest(args=args):
                result = self.run_cli("-", *args, input=">x\nA\n")
                self.assertEqual(result.returncode, 2)
                self.assertEqual(result.stdout, "")

    def test_gzip(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "input.FA.GZ"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(">x\nACGT\n")
            result = self.run_cli(path, "--format", "json")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads(result.stdout)["summary"]["length"], 4)

    def test_corrupt_gzip(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.gz"
            path.write_bytes(b"not gzip")
            result = self.run_cli(path)
            self.assertEqual(result.returncode, 2)
            self.assertNotIn("Traceback", result.stderr)

    def test_missing_file(self):
        result = self.run_cli("does-not-exist.fasta")
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, "")

    def test_invalid_utf8(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.fa"
            path.write_bytes(b">x\nAC\xff\n")
            self.assertEqual(self.run_cli(path).returncode, 2)

    def test_atomic_output_and_preservation_on_error(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "report.json"
            result = self.run_cli("-", "--format", "json", "-o", path, input=">x\nAC\n")
            self.assertEqual(result.returncode, 0, result.stderr)
            original = path.read_bytes()
            self.assertEqual(json.loads(original)["summary"]["length"], 2)
            result = self.run_cli("-", "-o", path, input=">x\nU\n")
            self.assertEqual(result.returncode, 2)
            self.assertEqual(path.read_bytes(), original)
            self.assertEqual(list(Path(folder).iterdir()), [path])

    def test_input_overwrite_blocked(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "input.fa"
            path.write_text(">x\nACGT\n", encoding="utf-8")
            original = path.read_bytes()
            result = self.run_cli(path, "-o", path)
            self.assertEqual(result.returncode, 2)
            self.assertEqual(path.read_bytes(), original)

    def test_hardlink_overwrite_blocked(self):
        with tempfile.TemporaryDirectory() as folder:
            source, target = Path(folder) / "source.fa", Path(folder) / "alias.fa"
            source.write_text(">x\nAC\n", encoding="utf-8")
            try:
                os.link(source, target)
            except OSError:
                self.skipTest("filesystem does not support hardlinks")
            self.assertEqual(self.run_cli(source, "-o", target).returncode, 2)
            self.assertEqual(source.read_text(encoding="utf-8"), ">x\nAC\n")

    def test_unwritable_output(self):
        with tempfile.TemporaryDirectory() as folder:
            result = self.run_cli("-", "-o", Path(folder) / "missing" / "report.md", input=">x\nA\n")
            self.assertEqual(result.returncode, 2)
            self.assertNotIn("Traceback", result.stderr)

    def test_version(self):
        result = self.run_cli("--version")
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), "fasta-lens 0.1.0")


if __name__ == "__main__":
    unittest.main()
