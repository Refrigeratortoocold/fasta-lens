"""Streaming DNA validation. No sequence data is sent to a service."""

from collections import Counter
from dataclasses import asdict, dataclass
import math
from typing import Iterable, Iterator

from . import __version__

IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
CANONICAL = frozenset("ACGT")


class FastaError(ValueError):
    """Invalid input, including a source line number when available."""


@dataclass(frozen=True)
class Record:
    identifier: str
    description: str
    length: int
    counts: dict[str, int]
    lowercase_bases: int


@dataclass(frozen=True)
class Thresholds:
    min_length: int = 1
    max_ambiguity_pct: float = 5.0

    def __post_init__(self):
        if isinstance(self.min_length, bool) or not isinstance(self.min_length, int) or self.min_length < 1:
            raise ValueError("min_length must be a positive integer")
        if not math.isfinite(self.max_ambiguity_pct) or not 0 <= self.max_ambiguity_pct <= 100:
            raise ValueError("max_ambiguity_pct must be between 0 and 100")


def parse_fasta(lines: Iterable[str]) -> Iterator[Record]:
    """Read strict, unaligned IUPAC DNA FASTA, retaining counts, not sequences.

    Blank lines are ignored. Whitespace within sequence lines, gaps, RNA U,
    comments and non-ASCII letters are rejected. Identifiers must be unique.
    """
    seen: set[str] = set()
    identifier = None
    description = ""
    counts: Counter[str] = Counter()
    lowercase = 0
    header_line = 0
    for number, raw in enumerate(lines, 1):
        line = raw.rstrip("\r\n")
        if number == 1:
            line = line.removeprefix("\ufeff")
        if not line.strip():
            continue
        if line.startswith(">"):
            if identifier is not None:
                if not counts:
                    raise FastaError(f"line {header_line}: empty sequence for {identifier!r}")
                yield Record(identifier, description, sum(counts.values()), dict(sorted(counts.items())), lowercase)
            header = line[1:].strip()
            if not header:
                raise FastaError(f"line {number}: empty FASTA header")
            if any(ord(char) < 32 or ord(char) == 127 for char in header):
                raise FastaError(f"line {number}: control character in header")
            identifier = header.split(maxsplit=1)[0]
            if identifier in seen:
                raise FastaError(f"line {number}: duplicate identifier {identifier!r}")
            seen.add(identifier)
            description = header[len(identifier):].strip()
            counts = Counter()
            lowercase = 0
            header_line = number
            continue
        if identifier is None:
            raise FastaError(f"line {number}: sequence before first FASTA header")
        for column, char in enumerate(line, 1):
            if not char.isascii() or char.upper() not in IUPAC_DNA:
                raise FastaError(f"line {number}, column {column}: invalid DNA character {char!r}")
            counts[char.upper()] += 1
            lowercase += int(char.islower())
    if identifier is None:
        raise FastaError("no FASTA records found")
    if not counts:
        raise FastaError(f"line {header_line}: empty sequence for {identifier!r}")
    yield Record(identifier, description, sum(counts.values()), dict(sorted(counts.items())), lowercase)


def n50(lengths: Iterable[int]) -> int:
    """Length at which descending records cover at least half of all bases."""
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    cumulative = 0
    for length in ordered:
        cumulative += length
        if cumulative * 2 >= total:
            return length
    return 0


def metrics(counts: dict[str, int]) -> dict:
    length = sum(counts.values())
    canonical = sum(counts.get(base, 0) for base in CANONICAL)
    ambiguous = length - canonical
    return {
        "length": length,
        "canonical_bases": canonical,
        "ambiguous_bases": ambiguous,
        "n_bases": counts.get("N", 0),
        "gc_pct": 100 * (counts.get("G", 0) + counts.get("C", 0)) / canonical if canonical else None,
        "ambiguity_pct": 100 * ambiguous / length if length else 0.0,
        "n_pct": 100 * counts.get("N", 0) / length if length else 0.0,
    }


def analyze(lines: Iterable[str], thresholds: Thresholds | None = None) -> dict:
    """Return a deterministic JSON-compatible report; raise on malformed input."""
    thresholds = thresholds or Thresholds()
    records = []
    totals: Counter[str] = Counter()
    for record in parse_fasta(lines):
        values = metrics(record.counts)
        warnings = []
        if record.length < thresholds.min_length:
            warnings.append("below_min_length")
        if values["ambiguity_pct"] > thresholds.max_ambiguity_pct:
            warnings.append("above_max_ambiguity")
        records.append({"id": record.identifier, "description": record.description,
                        **values, "lowercase_bases": record.lowercase_bases,
                        "base_counts": record.counts, "warnings": warnings})
        totals.update(record.counts)
    lengths = [record["length"] for record in records]
    warned = sum(bool(record["warnings"]) for record in records)
    return {
        "schema_version": 1,
        "tool": {"name": "fasta-lens", "version": __version__},
        "thresholds": asdict(thresholds),
        "summary": {"record_count": len(records), **metrics(dict(totals)),
                    "min_length": min(lengths), "max_length": max(lengths),
                    "mean_length": sum(lengths) / len(lengths), "n50": n50(lengths),
                    "warning_record_count": warned,
                    "lowercase_bases": sum(record["lowercase_bases"] for record in records)},
        "records": records,
    }
