from __future__ import annotations

import copy
import csv
import difflib
import io
import math
import re
import unicodedata
from collections import defaultdict

import numpy as np

DATE_COLUMN_HINTS = {"fecha", "fecha_adquisicion", "fecha_de_adquisicion", "date", "acquisition_date"}
MIN_CONTROL_N_RECOMMENDED = 3
SIMILAR_LABEL_THRESHOLD = 0.92


def normalize_label(label):
    text = "" if label is None else str(label)
    text = text.strip()
    text = unicodedata.normalize("NFKD", text)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = text.lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text


def tokenize_label(label):
    return {token for token in normalize_label(label).split("_") if token}


def uniquify_labels(labels):
    counts = {}
    result = {}
    for label in labels:
        base = normalize_label(label) or "unnamed_condition"
        counts[base] = counts.get(base, 0) + 1
        result[label] = base if counts[base] == 1 else f"{base}__{counts[base]}"
    return result


def parse_measurement(value):
    if value is None:
        return None, ""
    raw = str(value).strip()
    if raw == "":
        return None, raw
    lowered = raw.lower()
    if lowered in {"nan", "na", "n/a", "none", "null", "nd", "s/d"}:
        return None, raw
    clean = raw.replace("\u00a0", "").replace(" ", "")
    if "," in clean and "." in clean:
        if clean.rfind(",") > clean.rfind("."):
            clean = clean.replace(".", "").replace(",", ".")
        else:
            clean = clean.replace(",", "")
    elif "," in clean:
        clean = clean.replace(".", "").replace(",", ".")
    try:
        return float(clean), raw
    except ValueError:
        return None, raw


def parse_yy_mm_dd(raw_value):
    if raw_value is None:
        return None, ""
    raw = str(raw_value).strip()
    if raw == "":
        return None, raw
    if re.fullmatch(r"\d+\.0", raw):
        raw = raw[:-2]
    digits = re.sub(r"\D", "", raw)
    if len(digits) != 6:
        return None, raw
    from datetime import datetime
    try:
        parsed = datetime.strptime(digits, "%y%m%d").date()
        return parsed.isoformat(), raw
    except ValueError:
        return None, raw


def safe_mean(values):
    if not values:
        return None
    return float(np.mean(np.asarray(values, dtype=float)))


def safe_median(values):
    if not values:
        return None
    return float(np.median(np.asarray(values, dtype=float)))


def safe_std(values):
    if len(values) < 2:
        return None
    return float(np.std(np.asarray(values, dtype=float), ddof=1))


def safe_mad(values):
    if not values:
        return None
    arr = np.asarray(values, dtype=float)
    median = float(np.median(arr))
    return float(np.median(np.abs(arr - median)))


def safe_iqr(values):
    if not values:
        return None
    arr = np.asarray(values, dtype=float)
    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))
    return q3 - q1


def safe_cv(values):
    if len(values) < 2:
        return None
    mean = safe_mean(values)
    std = safe_std(values)
    if mean in (None, 0.0) or std is None:
        return None
    return std / mean


def compute_iqr_bounds(values):
    if len(values) < 4:
        return None
    arr = np.asarray(values, dtype=float)
    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))
    iqr = q3 - q1
    return q1 - 1.5 * iqr, q3 + 1.5 * iqr


def log2_or_none(value):
    if value is None or value <= 0:
        return None
    return math.log2(value)


def add_warning(warnings, level, code, message, **kwargs):
    warnings.append({
        "level": level,
        "code": code,
        "message": message,
        "source_file": kwargs.get("source_file", ""),
        "row_number": kwargs.get("row_number", ""),
        "column_name": kwargs.get("column_name", ""),
        "value": kwargs.get("value", ""),
        "date_raw": kwargs.get("date_raw", ""),
    })


def read_csv_text(text, source_file, warnings):
    sample = text[:4096]
    delimiter = ","
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
    except csv.Error:
        pass
    rows = list(csv.reader(io.StringIO(text), delimiter=delimiter))
    if not rows:
        raise ValueError(f"{source_file} is empty.")
    headers = [cell.strip() for cell in rows[0]]
    normalized_headers = []
    for index, name in enumerate(headers, start=1):
        if name:
            normalized_headers.append(name)
        else:
            placeholder = f"unnamed_column_{index}"
            normalized_headers.append(placeholder)
            add_warning(warnings, "warning", "blank_header", f"Blank header replaced with {placeholder}.", source_file=source_file, column_name=placeholder)
    parsed_rows = []
    width = len(normalized_headers)
    for row_index, row in enumerate(rows[1:], start=2):
        if not any(str(cell).strip() for cell in row):
            continue
        if len(row) < width:
            add_warning(warnings, "warning", "short_row", f"Row padded to {width} columns.", source_file=source_file, row_number=row_index)
            row = row + [""] * (width - len(row))
        elif len(row) > width:
            add_warning(warnings, "warning", "wide_row", "Extra cells were ignored.", source_file=source_file, row_number=row_index)
            row = row[:width]
        record = dict(zip(normalized_headers, row))
        record["__row_number__"] = str(row_index)
        parsed_rows.append(record)
    return {"headers": normalized_headers, "rows": parsed_rows}


def detect_date_column(headers):
    exact_matches, fuzzy_matches = [], []
    for header in headers:
        normalized = normalize_label(header)
        tokens = tokenize_label(header)
        if normalized in DATE_COLUMN_HINTS:
            exact_matches.append(header)
        elif "fecha" in tokens or "date" in tokens:
            fuzzy_matches.append(header)
    matches = exact_matches or fuzzy_matches
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f"Multiple date-like columns detected: {matches}")
    raise ValueError("No date column detected.")


def detect_default_control_column(headers, date_column):
    exact_matches, fuzzy_matches = [], []
    for header in headers:
        if header == date_column:
            continue
        normalized = normalize_label(header)
        tokens = tokenize_label(header)
        if normalized == "dmso":
            exact_matches.append(header)
        elif "dmso" in tokens or "dmso" in normalized:
            fuzzy_matches.append(header)
    matches = exact_matches or fuzzy_matches
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f"Multiple DMSO-like columns detected: {matches}")
    condition_headers = [header for header in headers if header != date_column]
    if condition_headers:
        return condition_headers[0]
    raise ValueError("No candidate control column detected.")


def warn_on_similar_columns(condition_headers, warnings, source_file):
    normalized_headers = {header: normalize_label(header) for header in condition_headers}
    by_normalized = defaultdict(list)
    for original, normalized in normalized_headers.items():
        by_normalized[normalized].append(original)
    for normalized, originals in by_normalized.items():
        if len(originals) > 1:
            add_warning(warnings, "warning", "duplicate_condition_normalized_name", f"Multiple columns collapse to '{normalized}': {originals}.", source_file=source_file)
    checked_pairs = set()
    for left in condition_headers:
        for right in condition_headers:
            if left == right:
                continue
            pair = tuple(sorted((left, right)))
            if pair in checked_pairs:
                continue
            checked_pairs.add(pair)
            ratio = difflib.SequenceMatcher(a=normalized_headers[left], b=normalized_headers[right]).ratio()
            if ratio >= SIMILAR_LABEL_THRESHOLD and normalized_headers[left] != normalized_headers[right]:
                add_warning(warnings, "warning", "similar_condition_names", f"Condition names may be inconsistent: {left!r} vs {right!r}.", source_file=source_file)


def build_long_rows(loaded, source_file, warnings, requested_control_column=None):
    headers = loaded["headers"]
    date_column = detect_date_column(headers)
    condition_headers = [header for header in headers if header != date_column]
    if requested_control_column:
        if requested_control_column not in condition_headers:
            raise ValueError(f"Selected control column {requested_control_column!r} is not present in {source_file}.")
        control_column = requested_control_column
    else:
        control_column = detect_default_control_column(headers, date_column)
    add_warning(warnings, "info", "detected_date_column", f"Detected date column: {date_column}.", source_file=source_file, column_name=date_column)
    add_warning(warnings, "info", "selected_control_column", f"Selected control column: {control_column}.", source_file=source_file, column_name=control_column)
    warn_on_similar_columns(condition_headers, warnings, source_file)
    unique_condition_keys = uniquify_labels(condition_headers)
    nonempty_counts = {header: 0 for header in condition_headers}
    long_rows = []
    for raw_row in loaded["rows"]:
        row_number = int(raw_row["__row_number__"])
        parsed_date, raw_date = parse_yy_mm_dd(raw_row.get(date_column))
        if parsed_date is None:
            add_warning(warnings, "warning", "invalid_date", "Row discarded because the date could not be parsed as YYMMDD.", source_file=source_file, row_number=row_number, column_name=date_column, value=raw_row.get(date_column, ""), date_raw=raw_row.get(date_column, ""))
            continue
        for condition_original in condition_headers:
            raw_value = raw_row.get(condition_original, "")
            numeric_value, cleaned_raw_value = parse_measurement(raw_value)
            if cleaned_raw_value != "":
                nonempty_counts[condition_original] += 1
            if cleaned_raw_value == "":
                continue
            if numeric_value is None:
                add_warning(warnings, "warning", "non_numeric_measurement", "Non-empty measurement discarded because it could not be parsed as numeric.", source_file=source_file, row_number=row_number, column_name=condition_original, value=cleaned_raw_value, date_raw=raw_date)
                continue
            long_rows.append({
                "source_file": source_file,
                "source_row_number": row_number,
                "date_raw": raw_date,
                "acquisition_date": parsed_date,
                "condition_original": condition_original,
                "condition_clean": normalize_label(condition_original),
                "condition_key": unique_condition_keys[condition_original],
                "raw_value": cleaned_raw_value,
                "intensity": float(numeric_value),
                "is_control_condition": condition_original == control_column,
            })
    for header, count in nonempty_counts.items():
        if count == 0:
            add_warning(warnings, "warning", "empty_condition_column", f"Condition column {header!r} produced no embryo records.", source_file=source_file, column_name=header)
    if not long_rows:
        raise ValueError("No valid embryo measurements were parsed from the file.")
    return long_rows, control_column


def summarize_values(values):
    return {
        "n_embryos": len(values),
        "mean": safe_mean(values),
        "median": safe_median(values),
        "sd": safe_std(values),
        "mad": safe_mad(values),
        "iqr": safe_iqr(values),
    }


def compute_outlier_flags(long_rows):
    groups = defaultdict(list)
    for row in long_rows:
        groups[(row["acquisition_date"], row["condition_key"])].append(row["intensity"])
    bounds_by_group = {key: compute_iqr_bounds(values) for key, values in groups.items()}
    for row in long_rows:
        bounds = bounds_by_group[(row["acquisition_date"], row["condition_key"])]
        if bounds is None:
            row["is_iqr_outlier_within_date_condition"] = False
            row["outlier_lower_bound"] = None
            row["outlier_upper_bound"] = None
            continue
        lower, upper = bounds
        row["outlier_lower_bound"] = lower
        row["outlier_upper_bound"] = upper
        row["is_iqr_outlier_within_date_condition"] = bool(row["intensity"] < lower or row["intensity"] > upper)


def normalize_by_control(long_rows, control_column, warnings, source_file, branch_label):
    control_by_date = defaultdict(list)
    for row in long_rows:
        if row["condition_original"] == control_column:
            control_by_date[row["acquisition_date"]].append(row["intensity"])
    control_rows = []
    control_summary = {}
    all_dates = sorted({row["acquisition_date"] for row in long_rows})
    for acquisition_date in all_dates:
        control_values = control_by_date.get(acquisition_date, [])
        stats = summarize_values(control_values)
        anchor = stats["median"]
        status = "ok"
        if not control_values:
            status = "missing_control"
            add_warning(warnings, "warning", "missing_control_anchor", f"Date has no valid control measurements in branch '{branch_label}'.", source_file=source_file, date_raw=acquisition_date, column_name=control_column)
        elif len(control_values) < MIN_CONTROL_N_RECOMMENDED:
            status = "low_control_n"
            add_warning(warnings, "warning", "low_control_n", f"Date has only {len(control_values)} control embryos in branch '{branch_label}'.", source_file=source_file, date_raw=acquisition_date, column_name=control_column)
        if anchor is None or anchor <= 0:
            status = "invalid_control_anchor"
            add_warning(warnings, "warning", "invalid_control_anchor", f"Date could not be normalized in branch '{branch_label}'.", source_file=source_file, date_raw=acquisition_date, column_name=control_column)
        control_row = {
            "source_file": source_file,
            "acquisition_date": acquisition_date,
            "control_condition_original": control_column,
            "control_n": len(control_values),
            "control_mean": stats["mean"],
            "control_median": anchor,
            "control_sd": stats["sd"],
            "control_mad": stats["mad"],
            "control_iqr": stats["iqr"],
            "anchor_status": status,
            "analysis_variant": branch_label,
        }
        control_rows.append(control_row)
        control_summary[acquisition_date] = control_row
    for row in long_rows:
        anchor_info = control_summary[row["acquisition_date"]]
        anchor = anchor_info["control_median"]
        row["control_n"] = anchor_info["control_n"]
        row["control_median"] = anchor
        row["anchor_status"] = anchor_info["anchor_status"]
        row["analysis_variant"] = branch_label
        if anchor is None or anchor <= 0:
            row["ratio_vs_control"] = None
            row["log2fc_vs_control"] = None
            row["normalization_status"] = "not_normalized_missing_anchor"
        else:
            ratio = row["intensity"] / anchor
            row["ratio_vs_control"] = ratio
            row["log2fc_vs_control"] = log2_or_none(ratio)
            row["normalization_status"] = "normalized"
    return long_rows, control_rows


def build_summary_rows(long_rows):
    grouped = defaultdict(list)
    for row in long_rows:
        grouped[(row["acquisition_date"], row["condition_key"])].append(row)
    summary_rows = []
    for (acquisition_date, condition_key), rows in sorted(grouped.items()):
        intensities = [row["intensity"] for row in rows]
        ratios = [row["ratio_vs_control"] for row in rows if row["ratio_vs_control"] is not None]
        log2fcs = [row["log2fc_vs_control"] for row in rows if row["log2fc_vs_control"] is not None]
        raw_stats = summarize_values(intensities)
        ratio_stats = summarize_values(ratios)
        log2_stats = summarize_values(log2fcs)
        summary_rows.append({
            "source_file": rows[0]["source_file"],
            "acquisition_date": acquisition_date,
            "condition_original": rows[0]["condition_original"],
            "condition_clean": rows[0]["condition_clean"],
            "condition_key": condition_key,
            "analysis_variant": rows[0]["analysis_variant"],
            "n_embryos": raw_stats["n_embryos"],
            "n_normalized": len(ratios),
            "n_outlier_flagged": sum(1 for row in rows if row["is_iqr_outlier_within_date_condition"]),
            "raw_mean": raw_stats["mean"],
            "raw_median": raw_stats["median"],
            "raw_sd": raw_stats["sd"],
            "raw_mad": raw_stats["mad"],
            "raw_iqr": raw_stats["iqr"],
            "normalized_ratio_mean": ratio_stats["mean"],
            "normalized_ratio_median": ratio_stats["median"],
            "normalized_ratio_sd": ratio_stats["sd"],
            "normalized_log2fc_mean": log2_stats["mean"],
            "normalized_log2fc_median": log2_stats["median"],
            "normalized_log2fc_sd": log2_stats["sd"],
            "anchor_status": rows[0]["anchor_status"],
        })
    return summary_rows


def build_variation_rows(summary_rows):
    by_condition = defaultdict(list)
    for row in summary_rows:
        by_condition[row["condition_key"]].append(row)
    variation_rows = []
    for condition_key, rows in sorted(by_condition.items()):
        raw_daily_medians = [row["raw_median"] for row in rows if row["raw_median"] is not None]
        normalized_daily_medians = [row["normalized_ratio_median"] for row in rows if row["normalized_ratio_median"] is not None]
        variation_rows.append({
            "source_file": rows[0]["source_file"],
            "condition_original": rows[0]["condition_original"],
            "condition_clean": rows[0]["condition_clean"],
            "condition_key": condition_key,
            "analysis_variant": rows[0]["analysis_variant"],
            "raw_daily_median_cv": safe_cv(raw_daily_medians),
            "normalized_daily_median_cv": safe_cv(normalized_daily_medians),
        })
    return variation_rows


def build_analysis_branch(base_rows, control_column, warnings, source_file, branch_label, exclude_outliers):
    branch_rows = [copy.deepcopy(row) for row in base_rows]
    if exclude_outliers:
        branch_rows = [row for row in branch_rows if not row["is_iqr_outlier_within_date_condition"]]
    if not branch_rows:
        add_warning(warnings, "warning", "empty_branch_after_outlier_removal", f"Branch '{branch_label}' has no rows after outlier exclusion.", source_file=source_file)
        return [], [], [], []
    branch_rows, control_rows = normalize_by_control(branch_rows, control_column, warnings, source_file, branch_label)
    branch_rows.sort(key=lambda row: (row["acquisition_date"], row["source_row_number"], row["condition_original"]))
    summary_rows = build_summary_rows(branch_rows)
    variation_rows = build_variation_rows(summary_rows)
    return branch_rows, summary_rows, control_rows, variation_rows


def rows_to_csv(rows):
    if not rows:
        return ""
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def prism_rows_to_csv(rows, headers):
    if not rows or not headers:
        return ""
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=headers, extrasaction="ignore", restval="")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _ordered_condition_keys(long_rows):
    """Returns condition keys in order: control first, then others in order of first appearance."""
    seen_keys = []
    seen_set = set()
    control_key = None
    key_to_original = {}
    for row in long_rows:
        ck = row["condition_key"]
        key_to_original[ck] = row["condition_original"]
        if ck not in seen_set:
            seen_set.add(ck)
            seen_keys.append(ck)
        if row["is_control_condition"]:
            control_key = ck
    ordered = []
    if control_key and control_key in seen_keys:
        ordered.append(control_key)
    for ck in seen_keys:
        if ck not in ordered:
            ordered.append(ck)
    return ordered, key_to_original


def build_prism_by_date(long_rows):
    """Wide Prism table: one row per acquisition date, one column per condition.
    Each cell = median log2FC of all embryos in that date-condition group.
    This is the recommended table for statistical inference (N = number of dates)."""
    if not long_rows:
        return [], []
    ordered_keys, key_to_original = _ordered_condition_keys(long_rows)
    groups = defaultdict(lambda: defaultdict(list))
    for row in long_rows:
        if row.get("log2fc_vs_control") is not None:
            groups[row["acquisition_date"]][row["condition_key"]].append(row["log2fc_vs_control"])
    all_dates = sorted({row["acquisition_date"] for row in long_rows})
    headers = ["acquisition_date"] + [key_to_original.get(ck, ck) for ck in ordered_keys]
    result_rows = []
    for date in all_dates:
        row_dict = {"acquisition_date": date}
        for ck in ordered_keys:
            col = key_to_original.get(ck, ck)
            values = groups[date].get(ck, [])
            row_dict[col] = safe_median(values) if values else ""
        result_rows.append(row_dict)
    return result_rows, headers


def build_prism_by_embryo(long_rows):
    """Wide Prism table: one column per condition, one row per individual embryo (all dates pooled).
    Each cell = log2FC of one embryo. Rows are NOT paired across conditions.
    Shorter columns are padded with empty cells.
    Use with caution: treating each embryo as independent may cause pseudo-replication."""
    if not long_rows:
        return [], []
    ordered_keys, key_to_original = _ordered_condition_keys(long_rows)
    by_condition = defaultdict(list)
    for row in long_rows:
        if row.get("log2fc_vs_control") is not None:
            by_condition[row["condition_key"]].append(row["log2fc_vs_control"])
    max_n = max((len(v) for v in by_condition.values()), default=0)
    if max_n == 0:
        return [], []
    headers = [key_to_original.get(ck, ck) for ck in ordered_keys]
    result_rows = []
    for i in range(max_n):
        row_dict = {}
        for ck in ordered_keys:
            col = key_to_original.get(ck, ck)
            values = by_condition.get(ck, [])
            row_dict[col] = values[i] if i < len(values) else ""
        result_rows.append(row_dict)
    return result_rows, headers


def analyze_one_file(source_file, text, control_column=None):
    warnings = []
    add_warning(warnings, "info", "file_loaded", "Archivo cargado en el navegador.", source_file=source_file)
    add_warning(warnings, "info", "raw_data_assumption",
        "SUPUESTO DE DATOS CRUDOS: Esta herramienta asume que todos los valores de entrada son "
        "intensidades de fluorescencia crudas sin procesar. Si tus datos ya fueron normalizados, "
        "escalados o transformados por el software de adquisición (p.ej. porcentaje de control, "
        "valores relativos), los resultados no serán válidos. Verifica el origen de tus datos "
        "antes de usar los resultados.",
        source_file=source_file)
    add_warning(warnings, "info", "multi_file_normalization_note",
        "NOTA MULTI-ARCHIVO: Cada archivo CSV es normalizado de forma completamente independiente, "
        "usando únicamente los embriones de control medidos en ese mismo archivo y en esa misma fecha. "
        "Si cargas varios archivos que comparten fechas de adquisición pero tienen embriones de control "
        "distintos (grupos diferentes), las condiciones entre archivos no son directamente comparables "
        "sin verificar que los anclas de normalización sean equivalentes. Interpreta las comparaciones "
        "entre archivos con precaución.",
        source_file=source_file)
    loaded = read_csv_text(text, source_file, warnings)
    long_rows, selected_control_column = build_long_rows(loaded, source_file, warnings, requested_control_column=control_column)
    compute_outlier_flags(long_rows)
    n_outliers = sum(1 for row in long_rows if row["is_iqr_outlier_within_date_condition"])
    if n_outliers:
        add_warning(warnings, "warning", "outliers_flagged_not_removed", f"Se detectaron {n_outliers} posibles outliers con la regla 1.5×IQR dentro de cada grupo fecha-condición.", source_file=source_file)
    add_warning(warnings, "info", "statistical_unit_note",
        "NOTA ESTADÍSTICA — PSEUDORREPLICACIÓN: La unidad experimental independiente es la FECHA de "
        "adquisición (cada día = una réplica experimental), no el embrión individual. "
        "La tabla 'PRISM_por_fecha' (una fila por día experimental) es la apropiada para inferencia "
        "estadística estándar (t-test, ANOVA, etc.) donde N = número de fechas. "
        "La tabla 'PRISM_por_embrion' trata cada embrión como observación independiente; "
        "usarla directamente en pruebas estadísticas puede producir pseudorreplicación y "
        "resultados con poder estadístico artificialmente inflado. El investigador decide qué "
        "tabla es apropiada para su diseño experimental.",
        source_file=source_file)
    removed_outlier_rows = [copy.deepcopy(row) for row in long_rows if row["is_iqr_outlier_within_date_condition"]]
    for row in removed_outlier_rows:
        row["analysis_variant"] = "removed_outliers"
    retained_long_rows, retained_summary_rows, retained_control_rows, retained_variation_rows = build_analysis_branch(long_rows, selected_control_column, warnings, source_file, "with_outliers", False)
    cleaned_long_rows, cleaned_summary_rows, cleaned_control_rows, cleaned_variation_rows = build_analysis_branch(long_rows, selected_control_column, warnings, source_file, "without_outliers", True)
    prism_date_with, prism_date_with_h = build_prism_by_date(retained_long_rows)
    prism_date_without, prism_date_without_h = build_prism_by_date(cleaned_long_rows)
    prism_embryo_with, prism_embryo_with_h = build_prism_by_embryo(retained_long_rows)
    prism_embryo_without, prism_embryo_without_h = build_prism_by_embryo(cleaned_long_rows)
    base = source_file.rsplit(".", 1)[0]
    output_files = {
        f"{base}_PRISM_por_fecha_sin_outliers.csv": prism_rows_to_csv(prism_date_without, prism_date_without_h),
        f"{base}_PRISM_por_fecha_con_outliers.csv": prism_rows_to_csv(prism_date_with, prism_date_with_h),
        f"{base}_PRISM_por_embrion_sin_outliers.csv": prism_rows_to_csv(prism_embryo_without, prism_embryo_without_h),
        f"{base}_PRISM_por_embrion_con_outliers.csv": prism_rows_to_csv(prism_embryo_with, prism_embryo_with_h),
        f"{base}_datos_normalizados_sin_outliers.csv": rows_to_csv(cleaned_long_rows),
        f"{base}_datos_normalizados_con_outliers.csv": rows_to_csv(retained_long_rows),
        f"{base}_outliers_eliminados.csv": rows_to_csv(removed_outlier_rows),
        f"{base}_resumen_por_fecha_condicion_sin_outliers.csv": rows_to_csv(cleaned_summary_rows),
        f"{base}_resumen_por_fecha_condicion_con_outliers.csv": rows_to_csv(retained_summary_rows),
        f"{base}_ancla_control_por_fecha_sin_outliers.csv": rows_to_csv(cleaned_control_rows),
        f"{base}_ancla_control_por_fecha_con_outliers.csv": rows_to_csv(retained_control_rows),
        f"{base}_variacion_por_condicion_sin_outliers.csv": rows_to_csv(cleaned_variation_rows),
        f"{base}_variacion_por_condicion_con_outliers.csv": rows_to_csv(retained_variation_rows),
        f"{base}_advertencias.csv": rows_to_csv(warnings),
    }
    return {
        "source_file": source_file,
        "summary": {
            "with_outliers_rows": len(retained_long_rows),
            "without_outliers_rows": len(cleaned_long_rows),
            "removed_outliers_rows": len(removed_outlier_rows),
            "warning_count": len(warnings),
            "n_dates": len(sorted({row["acquisition_date"] for row in retained_long_rows})),
            "n_conditions": len(sorted({row["condition_original"] for row in retained_long_rows})),
            "control_column": selected_control_column,
        },
        "warnings": warnings,
        "outputs": output_files,
        "branches": {
            "with_outliers": {"long_rows": retained_long_rows, "summary_rows": retained_summary_rows, "control_rows": retained_control_rows, "variation_rows": retained_variation_rows},
            "without_outliers": {"long_rows": cleaned_long_rows, "summary_rows": cleaned_summary_rows, "control_rows": cleaned_control_rows, "variation_rows": cleaned_variation_rows},
        },
    }


def analyze_file_payload(payload_json, progress_cb=None):
    import json
    payload = json.loads(payload_json)
    files = payload.get("files", [])
    control_overrides = payload.get("control_overrides", {})
    results = []
    ignored_files = payload.get("ignored_files", [])
    total = max(1, len(files))
    for index, item in enumerate(files):
        if progress_cb is not None:
            progress_cb(f"Analyzing {item['name']} ({index + 1}/{total})…", int(10 + 80 * index / total))
        results.append(analyze_one_file(item["name"], item["text"], control_overrides.get(item["name"])))
    if progress_cb is not None:
        progress_cb("Finalizing browser outputs…", 96)
    return {
        "results": results,
        "ignored_files": ignored_files,
        "run_summary": [{"source_file": item["source_file"], **item["summary"]} for item in results],
    }
