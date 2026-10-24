# FASTA Lens

[![Tests](https://github.com/Refrigeratortoocold/fasta-lens/actions/workflows/ci.yml/badge.svg)](https://github.com/Refrigeratortoocold/fasta-lens/actions/workflows/ci.yml)

**Check DNA FASTA files before using them in an analysis.**

FASTA Lens is a small, offline Python command line tool for students, educators,
and bioinformatics workflows. It validates unaligned DNA, finds ambiguous bases,
and produces reproducible reports without runtime dependencies or API keys.

**Status:** initial version 0.1.0. This is a new project maintained by Soh Jia Sheng.
External adoption and independent validation have not yet been established.
Development was assisted by OpenAI Codex.

## Quick start

Python 3.10 or newer is required.

```sh
git clone https://github.com/Refrigeratortoocold/fasta-lens.git
cd fasta-lens
python -m fasta_lens examples/demo.fasta
```

The module runs directly from this checkout with the Python standard library.
To install the command in a virtual environment:

```sh
python -m venv .venv
# macOS / Linux:
. .venv/bin/activate
# Windows PowerShell instead:
# .venv\Scripts\Activate.ps1
python -m pip install .
fasta-lens examples/demo.fasta
```

Installation uses setuptools as a build dependency. No PyPI release is claimed.

## Example

All included sequences are synthetic teaching examples.

| Record | Length (bp) | GC (%) | Ambiguous (%) | Warnings |
| --- | ---: | ---: | ---: | --- |
| clean | 12 | 50.00 | 0.00 | none |
| ambiguous | 8 | 50.00 | 50.00 | above_max_ambiguity |
| gc_rich | 10 | 80.00 | 0.00 | none |

Across these three records: **30 bases, N50 10 bp, GC 61.54%**.
See [the complete Markdown report](examples/demo-report.md) and
[JSON output](examples/demo-report.json).

## What it does

- Validates headers, unique IDs, nonempty sequences, and IUPAC DNA characters.
- Reads wrapped lines, lowercase bases, UTF-8 BOM, gzip files, and standard input.
- Reports length, canonical GC percentage, N percentage, ambiguity, and lowercase counts.
- Summarizes total bases, record counts, length range, mean length, and N50.
- Writes Markdown, JSON, or TSV with configurable per-record warning thresholds.
- Provides predictable exit codes for local checks and continuous integration.
- Keeps processing offline. There are no telemetry, network, or AI calls.

## Reports and thresholds

```sh
python -m fasta_lens reads.fasta.gz --format json -o report.json
python -m fasta_lens reads.fasta --format tsv -o report.tsv
python -m fasta_lens reads.fasta --min-length 100 --max-ambiguity-pct 2 --fail-on-warning
```

For piped input, use `-` as the input path. For example on macOS / Linux:

```sh
cat reads.fasta | python -m fasta_lens - --format json
```

| Exit code | Meaning |
| --- | --- |
| 0 | Valid input; warnings allowed unless requested otherwise |
| 1 | Valid input with warnings and `--fail-on-warning` enabled; report still written |
| 2 | Invalid input, invalid options, or a read/write error |

Default minimum length is 1 bp and maximum ambiguity is 5%.
Length equal to the minimum passes. Ambiguity equal to the maximum passes.
Warnings do not modify or remove input records. Invalid input produces no report.
Output files are replaced atomically after successful analysis; input-file
overwrites, including hardlink aliases, are blocked.

## Metric definitions

| Metric | Definition |
| --- | --- |
| Length | All valid IUPAC DNA bases, including ambiguous bases |
| GC (%) | 100 x (G + C) / (A + C + G + T) |
| Ambiguity (%) | 100 x non-ACGT bases / total bases; includes N |
| N (%) | 100 x N bases / total bases |
| Lowercase bases | Original lowercase symbols, including ambiguous symbols |
| N50 | Record length at which descending lengths cover at least half of all bases |

GC is `null` in JSON and `NA` in tabular output when a record has no A/C/G/T.
Summary percentages use pooled base counts, not averages of per-record percentages.
The accepted DNA alphabet is `ACGTRYSWKMBDHVN`, case insensitive.

N50 describes the lengths of the records supplied. It does not establish
assembly correctness or biological quality.

## Python API

```python
from fasta_lens.core import Thresholds, analyze

with open("reads.fasta", encoding="utf-8-sig") as handle:
    report = analyze(handle, Thresholds(min_length=100, max_ambiguity_pct=2))
print(report["summary"])
```

JSON includes `schema_version: 1`, the tool version, thresholds, summary, and
per-record results. Input order is preserved. Reports contain no timestamps or
absolute input paths, making repeated runs on identical input deterministic.

## Input policy and limits

This version is deliberately strict: it accepts unaligned DNA FASTA, not RNA,
proteins, FASTQ, or alignment files. It rejects U, gaps, digits, comments,
internal sequence whitespace, duplicate IDs, and control characters in headers.
Blank lines are ignored. A header begins with `>` in column one; its first
whitespace-delimited token is the record ID.

The parser reads one line at a time and retains base counts rather than full
sequences. Reports and the set of IDs remain in memory, so memory grows with
record count and the longest input line. Gzip is detected by a `.gz` extension.
TSV contains raw record IDs; treat untrusted IDs as text when importing into a
spreadsheet.

This is a basic data check, not an assessment of contamination, sequencing
accuracy, protein function, or clinical suitability. It does not filter or
repair sequences. The project is not affiliated with OpenAI or NUS.

## Development and contributions

```sh
python -m unittest discover -s tests -v
```

CI runs tests and an installed-command smoke check on Linux and Windows with
Python 3.10 and 3.13. See [CONTRIBUTING.md](CONTRIBUTING.md) for changes and
[ROADMAP.md](ROADMAP.md) for open work. Bug reports with small synthetic examples
are welcome.

Licensed under the [MIT License](LICENSE).
