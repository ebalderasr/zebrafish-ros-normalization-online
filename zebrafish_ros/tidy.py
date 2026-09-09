"""Reading the raw DCF fluorescence CSVs and converting them to tidy format.

Expected input layout (one wide matrix per file)::

    FECHA,DMSO,VAS,DPI,APO
    260321,118.4,102.7,96.2,99.1
    260321,124.9,98.3,,101.5
    ...

Each row is one embryo and each column one condition. Empty cells are embryos
that were not measured under that condition. Rows are not paired: two values on
the same row are two different embryos acquired on the same date, not one
embryo measured twice.

This is the same input shape that `hyper-normalizer` accepts, so a single
experiment can be exported once and processed by either pipeline. The file name
may follow the ``{GROUP} DCF {PANEL}.csv`` convention, in which case files that
share a group are analysed together and their common control is de-duplicated.
A file whose name does not match is treated as its own group.

The parser is deliberately permissive, matching the browser engine: it accepts
decimal commas, accented and inconsistently capitalised headers, numbers stored
as text, and blank columns.
"""

from __future__ import annotations

import csv
import re
import unicodedata
from dataclasses import dataclass
from collections import defaultdict
from datetime import datetime
from pathlib import Path

#: ``DCL DCF ANTINOX.csv`` -> group ``DCL``, panel ``ANTINOX``.
FILENAME_RE = re.compile(r"^(?P<group>[^ ]+) +\S+ +(?P<panel>[^ ]+)\.csv$", re.IGNORECASE)

#: Header names recognised as the acquisition date column.
DATE_COLUMN_HINTS = {"fecha", "fecha_adquisicion", "fecha_de_adquisicion", "date", "acquisition_date"}

#: Cell contents that mean "no measurement".
NULL_TOKENS = {"nan", "na", "n/a", "none", "null", "nd", "s/d"}

DEFAULT_CONTROL = "DMSO"


@dataclass(frozen=True)
class Measurement:
    """One observation: a single embryo under one condition on one date."""

    group: str          # experiment group, from the file name
    panel: str          # drug panel, from the file name; empty when absent
    date: str           # acquisition date as ISO yyyy-mm-dd
    date_raw: str       # the date exactly as written in the file
    treatment: str      # condition column header
    intensity: float    # raw DCF fluorescence intensity
    source: str         # originating file, kept for traceability
    row_number: int     # 1-based row in the source file


class InputError(ValueError):
    """The input directory or files do not have the expected shape."""


def normalize_label(label: str | None) -> str:
    """Fold a header to a comparable form: no accents, no case, no padding."""
    text = "" if label is None else str(label)
    text = unicodedata.normalize("NFKD", text)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    return re.sub(r"\s+", " ", text).strip().lower()


def parse_measurement(value: str | None) -> float | None:
    """Parse one cell, accepting decimal commas and thousands separators."""
    if value is None:
        return None
    raw = str(value).strip()
    if raw == "" or raw.lower() in NULL_TOKENS:
        return None

    clean = raw.replace(" ", "").replace(" ", "")
    if "," in clean and "." in clean:
        # Whichever separator comes last is the decimal one.
        if clean.rfind(",") > clean.rfind("."):
            clean = clean.replace(".", "").replace(",", ".")
        else:
            clean = clean.replace(",", "")
    elif "," in clean:
        clean = clean.replace(".", "").replace(",", ".")
    try:
        return float(clean)
    except ValueError:
        return None


def parse_date(value: str | None) -> str | None:
    """Parse a YYMMDD date into ISO form, or return ``None``.

    ISO dates sort chronologically as text, which YYMMDD does not once a
    century boundary is crossed.
    """
    if value is None:
        return None
    raw = str(value).strip()
    if raw == "":
        return None
    if re.fullmatch(r"\d+\.0", raw):  # a spreadsheet turned the date into a float
        raw = raw[:-2]
    digits = re.sub(r"\D", "", raw)
    if len(digits) != 6:
        return None
    try:
        return datetime.strptime(digits, "%y%m%d").date().isoformat()
    except ValueError:
        return None


def parse_filename(path: Path) -> tuple[str, str]:
    """Return ``(group, panel)`` for a file, falling back to its stem."""
    match = FILENAME_RE.match(path.name)
    if match is None:
        return path.stem, ""
    return match["group"].upper(), match["panel"].upper()


def discover_inputs(directory: Path) -> list[Path]:
    """Return every CSV in the directory, sorted."""
    paths = sorted(directory.glob("*.csv"))
    if not paths:
        raise InputError(f"No CSV file found in {directory}")
    return paths


def detect_date_column(headers: list[str]) -> str:
    """Find the acquisition date column by name, then by content shape."""
    for header in headers:
        if normalize_label(header).replace(" ", "_") in DATE_COLUMN_HINTS:
            return header
    if headers:
        return headers[0]
    raise InputError("The file has no header row.")


def read_wide(path: Path) -> list[Measurement]:
    """Read one wide CSV and return a measurement per non-empty numeric cell."""
    group, panel = parse_filename(path)

    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise InputError(f"{path.name} is empty") from None

        headers = [h.strip() for h in header]
        date_column = detect_date_column(headers)
        date_index = headers.index(date_column)

        out: list[Measurement] = []
        for row_number, row in enumerate(reader, start=2):
            if not row:
                continue
            date = parse_date(row[date_index] if date_index < len(row) else None)
            if date is None:
                continue
            date_raw = (row[date_index] if date_index < len(row) else "").strip()
            for index, treatment in enumerate(headers):
                if index == date_index:
                    continue
                intensity = parse_measurement(row[index] if index < len(row) else None)
                if intensity is None:
                    continue
                out.append(
                    Measurement(
                        group=group,
                        panel=panel,
                        date=date,
                        date_raw=date_raw,
                        treatment=treatment,
                        intensity=intensity,
                        source=path.name,
                        row_number=row_number,
                    )
                )
    return out


def deduplicate_controls(
    measurements: list[Measurement], control: str = DEFAULT_CONTROL
) -> tuple[list[Measurement], list[str]]:
    """Keep a single copy of the control per ``(group, date)``.

    When one experiment is split across drug panels, the control is measured
    once and copied into each panel's sheet. Counting it once per file would
    multiply its n and bias the normalization anchor. Files that do not share a
    group are unaffected.
    """
    kept = [m for m in measurements if m.treatment != control]
    warnings: list[str] = []

    by_group: dict[tuple[str, str], dict[str, list[Measurement]]] = defaultdict(dict)
    for m in measurements:
        if m.treatment == control:
            by_group[(m.group, m.date)].setdefault(m.source, []).append(m)

    for (group, date), per_source in sorted(by_group.items()):
        # Most replicates wins. Ties go to the alphabetically first file, so the
        # attribution in the SOURCE column does not depend on directory order.
        winner_source, winner = min(
            per_source.items(), key=lambda item: (-len(item[1]), item[0])
        )
        winner_values = [m.intensity for m in winner]

        for source, series in sorted(per_source.items()):
            if source == winner_source:
                continue
            values = [m.intensity for m in series]
            shared = min(len(values), len(winner_values))
            if values[:shared] != winner_values[:shared]:
                warnings.append(
                    f"Control '{control}' for {group}/{date} differs between "
                    f"'{winner_source}' and '{source}'. Used '{winner_source}'. "
                    "Check the data entry: both should hold the same control."
                )
        kept.extend(winner)

    return kept, warnings


def build_tidy(
    paths: list[Path], control: str = DEFAULT_CONTROL
) -> tuple[list[Measurement], list[str]]:
    """Read every file and return the tidy dataset with a single control series."""
    raw: list[Measurement] = []
    for path in paths:
        raw.extend(read_wide(path))
    if not raw:
        raise InputError("No valid embryo measurements were parsed.")
    return deduplicate_controls(raw, control=control)
