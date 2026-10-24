import io
import unittest
from fasta_lens.core import FastaError, Thresholds, analyze, n50, parse_fasta


class CoreTests(unittest.TestCase):
    def report(self, text, thresholds=None):
        return analyze(io.StringIO(text), thresholds)

    def test_wrapping_case_and_counts(self):
        record = self.report(">one description\nAcgT\nGC\n")["records"][0]
        self.assertEqual(record["length"], 6)
        self.assertEqual(record["lowercase_bases"], 2)
        self.assertEqual(record["description"], "description")
        self.assertEqual(record["base_counts"], {"A": 1, "C": 2, "G": 2, "T": 1})
        self.assertAlmostEqual(record["gc_pct"], 100 * 4 / 6)

    def test_all_iupac_symbols(self):
        record = self.report(">x\nACGTRYSWKMBDHVN\n")["records"][0]
        self.assertEqual(record["ambiguous_bases"], 11)
        self.assertEqual(record["n_bases"], 1)
        self.assertEqual(record["gc_pct"], 50)

    def test_all_ambiguous_has_no_gc(self):
        summary = self.report(">x\nNNRY\n")["summary"]
        self.assertIsNone(summary["gc_pct"])
        self.assertEqual(summary["ambiguity_pct"], 100)
        self.assertEqual(summary["n_pct"], 50)

    def test_weighted_gc(self):
        summary = self.report(">a\nG\n>b\nAAAAAAAAAN\n")["summary"]
        self.assertEqual(summary["gc_pct"], 10)
        self.assertEqual(summary["length"], 11)

    def test_bom_blank_lines_crlf_and_no_final_newline(self):
        report = self.report("\ufeff>a\r\nAC\r\n\r\nGT")
        self.assertEqual(report["summary"]["length"], 4)

    def test_invalid_inputs(self):
        cases = [
            ("", "no FASTA records"),
            (" \n\n", "no FASTA records"),
            ("ACGT\n", "line 1"),
            (">\nAC\n", "empty FASTA header"),
            (">a\n", "empty sequence"),
            (">a\n>b\nA\n", "empty sequence"),
            (">a\nA\n>a other\nT\n", "duplicate identifier"),
            (">a\nAC U\n", "column 3"),
            (">a\nACU\n", "column 3"),
            (">a\nAC-\n", "column 3"),
            (">a\nAC1\n", "column 3"),
            (">a\nA?\n", "column 2"),
            (">a\tname\nAC\n", "control character"),
            (">a\n;comment\n", "invalid DNA character"),
            (" >a\nAC\n", "sequence before first"),
        ]
        for text, message in cases:
            with self.subTest(text=text), self.assertRaisesRegex(FastaError, message):
                self.report(text)

    def test_threshold_boundaries(self):
        report = self.report(">a\nACGN\n", Thresholds(4, 25))
        self.assertEqual(report["records"][0]["warnings"], [])
        report = self.report(">a\nACGN\n", Thresholds(5, 24))
        self.assertEqual(report["records"][0]["warnings"], ["below_min_length", "above_max_ambiguity"])
        self.assertEqual(report["summary"]["warning_record_count"], 1)

    def test_bad_thresholds(self):
        for value in (-1, 101, float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                Thresholds(max_ambiguity_pct=value)
        for value in (0, -1, 1.5, True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                Thresholds(min_length=value)

    def test_n50(self):
        for lengths, expected in [([], 0), ([5], 5), ([8, 4, 4], 8), ([6, 5, 4], 5), ([1, 1, 1], 1)]:
            with self.subTest(lengths=lengths):
                self.assertEqual(n50(lengths), expected)

    def test_multi_record_summary(self):
        summary = self.report(">a\nAAAAAA\n>b\nGGGGG\n>c\nTTTT\n")["summary"]
        self.assertEqual((summary["record_count"], summary["min_length"], summary["max_length"]), (3, 4, 6))
        self.assertEqual(summary["n50"], 5)
        self.assertEqual(summary["mean_length"], 5)

    def test_deterministic(self):
        text = ">a\nACGTRYN\n"
        self.assertEqual(self.report(text), self.report(text))

    def test_parser_yields_before_reading_next_sequence(self):
        lines = iter([">a\n", "AC\n", ">b\n", "BAD!\n"])
        records = parse_fasta(lines)
        self.assertEqual(next(records).identifier, "a")
        with self.assertRaises(FastaError):
            next(records)

    def test_single_pass_iterable(self):
        self.assertEqual(analyze(iter([">a\n", "ACGT\n"]))["summary"]["length"], 4)


if __name__ == "__main__":
    unittest.main()
