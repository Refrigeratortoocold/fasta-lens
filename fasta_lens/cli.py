"""Command line interface and portable report renderers."""

import argparse
import csv
import gzip
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import zlib

from . import __version__
from .core import FastaError, Thresholds, analyze


def _number(value):
    return "NA" if value is None else f"{value:.2f}"


def _markdown_cell(value):
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\\", "&#92;").replace("|", "&#124;").replace("`", "&#96;").replace("*", "&#42;").replace("_", "&#95;").replace("[", "&#91;").replace("]", "&#93;")


def render(report: dict, format_name: str) -> str:
    if format_name == "json":
        return json.dumps(report, indent=2, ensure_ascii=True, allow_nan=False) + "\n"
    if format_name == "tsv":
        output = io.StringIO(newline="")
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "length", "gc_pct", "ambiguity_pct", "n_pct", "lowercase_bases", "warnings"])
        for record in report["records"]:
            writer.writerow([record["id"], record["length"], _number(record["gc_pct"]),
                             _number(record["ambiguity_pct"]), _number(record["n_pct"]),
                             record["lowercase_bases"], ",".join(record["warnings"])])
        return output.getvalue()
    if format_name != "markdown":
        raise ValueError(f"unknown format: {format_name}")
    summary = report["summary"]
    lines = ["# FASTA Lens quality report", "",
             f"Records: {summary['record_count']} | Bases: {summary['length']} | N50: {summary['n50']} bp",
             f"GC: {_number(summary['gc_pct'])}% | Ambiguous bases: {_number(summary['ambiguity_pct'])}%",
             f"Records with warnings: {summary['warning_record_count']}", "",
             "| Record | Length (bp) | GC (%) | Ambiguous (%) | Warnings |",
             "| --- | ---: | ---: | ---: | --- |"]
    for record in report["records"]:
        lines.append(f"| {_markdown_cell(record['id'])} | {record['length']} | {_number(record['gc_pct'])} | {_number(record['ambiguity_pct'])} | {', '.join(record['warnings']) or 'none'} |")
    lines += ["", "GC uses A+C+G+T as its denominator. NA means no canonical bases.",
              "Ambiguity includes all non-ACGT IUPAC DNA symbols, including N.", ""]
    return "\n".join(lines)


def _write_atomic(path: Path, content: str):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Validate unaligned DNA FASTA and report quality metrics offline.")
    parser.add_argument("input", help="UTF-8 FASTA file, .gz file, or - for stdin")
    parser.add_argument("--format", choices=("json", "tsv", "markdown"), default="markdown")
    parser.add_argument("-o", "--output", type=Path, help="write report atomically to this file")
    parser.add_argument("--min-length", type=int, default=1, metavar="BP")
    parser.add_argument("--max-ambiguity-pct", type=float, default=5.0, metavar="PCT")
    parser.add_argument("--fail-on-warning", action="store_true", help="exit 1 if any record exceeds a threshold")
    parser.add_argument("--version", action="version", version=f"fasta-lens {__version__}")
    args = parser.parse_args(argv)
    try:
        thresholds = Thresholds(args.min_length, args.max_ambiguity_pct)
        if args.output and args.input != "-":
            source, target = Path(args.input), args.output
            if source.resolve() == target.resolve() or (target.exists() and os.path.samefile(source, target)):
                raise ValueError("output must not overwrite the input file")
        if args.input == "-":
            report = analyze(sys.stdin, thresholds)
        else:
            opener = gzip.open if args.input.lower().endswith(".gz") else open
            with opener(args.input, "rt", encoding="utf-8-sig") as handle:
                report = analyze(handle, thresholds)
        content = render(report, args.format)
        if args.output:
            _write_atomic(args.output, content)
        else:
            sys.stdout.write(content)
        return int(args.fail_on_warning and report["summary"]["warning_record_count"] > 0)
    except (FastaError, OSError, UnicodeError, ValueError, EOFError, zlib.error) as error:
        print(f"fasta-lens: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
