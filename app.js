"use strict";

let pyodide = null;
let selectedFiles = [];
let ignoredFiles = [];
let analysisResults = null;

let elInitSection, elInitBar, elInitMsg;
let elUploadSection, elDropZone, elFileInput, elAnalyzeBtn, elAnalyzeHint;
let elProcessSection, elProcessBar, elProcessMsg;
let elSelectedFilesList, elSelectedFilesCount;
let elIgnoredFilesBox, elIgnoredFilesList;
let elResultsSection, elRunSummaryCards, elFileTabs, elFileTabContent, elDownloadAllBtn;

const SET2 = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"];
const MARKERS = ["circle", "square", "triangle-up", "diamond", "cross", "x", "triangle-down", "triangle-left", "triangle-right", "star"];

document.addEventListener("DOMContentLoaded", () => {
  elInitSection = document.getElementById("init-section");
  elInitBar = document.getElementById("init-bar");
  elInitMsg = document.getElementById("init-msg");
  elUploadSection = document.getElementById("upload-section");
  elDropZone = document.getElementById("drop-zone");
  elFileInput = document.getElementById("file-input");
  elAnalyzeBtn = document.getElementById("analyze-btn");
  elAnalyzeHint = document.getElementById("analyze-hint");
  elProcessSection = document.getElementById("process-section");
  elProcessBar = document.getElementById("process-bar");
  elProcessMsg = document.getElementById("process-msg");
  elSelectedFilesList = document.getElementById("selected-files-list");
  elSelectedFilesCount = document.getElementById("selected-files-count");
  elIgnoredFilesBox = document.getElementById("ignored-files-box");
  elIgnoredFilesList = document.getElementById("ignored-files-list");
  elResultsSection = document.getElementById("results-section");
  elRunSummaryCards = document.getElementById("run-summary-cards");
  elFileTabs = document.getElementById("file-tabs");
  elFileTabContent = document.getElementById("file-tab-content");
  elDownloadAllBtn = document.getElementById("download-all-btn");

  setupDropZone();
  elAnalyzeBtn.addEventListener("click", analyzeSelectedFiles);
  elDownloadAllBtn.addEventListener("click", downloadAllOutputsZip);
  initPyodide();
});

async function initPyodide() {
  try {
    setInitProgress("Loading Python runtime…", 15);
    setInitStep(1, "active");
    pyodide = await loadPyodide();
    setInitProgress("Installing scientific packages…", 42);
    setInitStep(1, "done");
    setInitStep(2, "active");
    await pyodide.loadPackage(["numpy"]);
    setInitProgress("Loading zebrafish ROS engine…", 78);
    setInitStep(2, "done");
    setInitStep(3, "active");
    const resp = await fetch("zebrafish_ros_engine.py");
    if (!resp.ok) throw new Error("Could not load zebrafish_ros_engine.py");
    const code = await resp.text();
    await pyodide.runPythonAsync(code);
    setInitProgress("Ready!", 100);
    setInitStep(3, "done");
    await sleep(450);
    show(elUploadSection);
    hide(elInitSection);
  } catch (err) {
    elInitMsg.textContent = "Error loading Python environment: " + err.message;
    elInitMsg.style.color = "var(--md-error)";
    console.error(err);
  }
}

function setInitStep(n, state) {
  const step = document.getElementById(`init-step-${n}`);
  if (!step) return;
  step.classList.remove("active", "done");
  step.classList.add(state);
  const dot = step.querySelector(".init-step-dot");
  if (dot && state === "done") dot.textContent = "✓";
}

function setInitProgress(msg, pct) {
  elInitMsg.textContent = msg;
  elInitBar.style.width = `${pct}%`;
}

function setProcessProgress(msg, pct) {
  elProcessMsg.textContent = msg;
  elProcessBar.style.width = `${pct}%`;
}

function setupDropZone() {
  elDropZone.addEventListener("dragover", (e) => {
    e.preventDefault();
    elDropZone.classList.add("drag-over");
  });
  elDropZone.addEventListener("dragleave", () => {
    elDropZone.classList.remove("drag-over");
  });
  elDropZone.addEventListener("drop", (e) => {
    e.preventDefault();
    elDropZone.classList.remove("drag-over");
    handleIncomingFiles(Array.from(e.dataTransfer.files || []));
  });
  elDropZone.addEventListener("click", () => elFileInput.click());
  elFileInput.addEventListener("change", (e) => {
    handleIncomingFiles(Array.from(e.target.files || []));
  });
}

function handleIncomingFiles(files) {
  const csvFiles = [];
  const ignored = [];
  for (const file of files) {
    if (file.name.toLowerCase().endsWith(".csv")) csvFiles.push(file);
    else ignored.push(file.name);
  }
  selectedFiles = csvFiles;
  ignoredFiles = ignored;
  renderSelectedFiles();
}

function renderSelectedFiles() {
  elSelectedFilesCount.textContent = `${selectedFiles.length} CSV file${selectedFiles.length === 1 ? "" : "s"}`;
  if (!selectedFiles.length) {
    elSelectedFilesList.classList.add("empty");
    elSelectedFilesList.textContent = "No files selected yet.";
    elAnalyzeBtn.disabled = true;
    elAnalyzeHint.textContent = "Select one or more CSV files to continue";
  } else {
    elSelectedFilesList.classList.remove("empty");
    elSelectedFilesList.innerHTML = selectedFiles.map((file) => `
      <div class="selected-file-item">
        <span class="file-chip-name">${escapeHtml(file.name)}</span>
        <span class="file-chip-type">${formatBytes(file.size)}</span>
      </div>
    `).join("");
    elAnalyzeBtn.disabled = false;
    elAnalyzeHint.textContent = "Ready to analyze in the browser";
  }
  if (ignoredFiles.length) {
    elIgnoredFilesBox.classList.remove("d-none");
    elIgnoredFilesList.innerHTML = ignoredFiles.map((name) => `<div class="ignored-item">${escapeHtml(name)}</div>`).join("");
  } else {
    elIgnoredFilesBox.classList.add("d-none");
    elIgnoredFilesList.innerHTML = "";
  }
}

async function analyzeSelectedFiles() {
  if (!selectedFiles.length) return;
  hide(elUploadSection);
  hide(elResultsSection);
  show(elProcessSection);
  setProcessProgress("Reading selected files…", 5);
  try {
    const files = [];
    for (const file of selectedFiles) {
      files.push({ name: file.name, text: await file.text() });
    }
    pyodide.globals.set("_progress_cb", (msg, pct) => setProcessProgress(msg, pct));
    pyodide.globals.set("_payload_json", JSON.stringify({ files, ignored_files: ignoredFiles }));
    const result = await pyodide.runPythonAsync("analyze_file_payload(_payload_json, _progress_cb)");
    analysisResults = deepConvert(result);
    setProcessProgress("Rendering results…", 100);
    renderResults(analysisResults);
    hide(elProcessSection);
    show(elResultsSection);
    elResultsSection.scrollIntoView({ behavior: "smooth" });
  } catch (err) {
    setProcessProgress("Error: " + err.message, 0);
    console.error(err);
    setTimeout(() => {
      hide(elProcessSection);
      show(elUploadSection);
    }, 2500);
  }
}

function renderResults(payload) {
  renderRunSummary(payload.run_summary || []);
  renderFileTabs(payload.results || []);
}

function renderRunSummary(summaryRows) {
  const totalFiles = summaryRows.length;
  const totalWarnings = summaryRows.reduce((sum, row) => sum + Number(row.warning_count || 0), 0);
  const totalRemoved = summaryRows.reduce((sum, row) => sum + Number(row.removed_outliers_rows || 0), 0);
  const totalWith = summaryRows.reduce((sum, row) => sum + Number(row.with_outliers_rows || 0), 0);
  const totalWithout = summaryRows.reduce((sum, row) => sum + Number(row.without_outliers_rows || 0), 0);
  elRunSummaryCards.innerHTML = [
    summaryCard("CSV files", totalFiles),
    summaryCard("Rows with outliers", totalWith),
    summaryCard("Rows without outliers", totalWithout),
    summaryCard("Removed outliers", totalRemoved),
    summaryCard("Warnings", totalWarnings),
  ].join("");
}

function summaryCard(label, value) {
  return `<div class="summary-stat-card"><div class="summary-stat-label">${escapeHtml(label)}</div><div class="summary-stat-value">${escapeHtml(String(value))}</div></div>`;
}

function renderFileTabs(results) {
  elFileTabs.innerHTML = results.map((item, index) => `
    <li class="nav-item" role="presentation">
      <button class="nav-link ${index === 0 ? "active" : ""}" id="file-tab-${index}" data-bs-toggle="tab" data-bs-target="#file-pane-${index}" type="button" role="tab">
        ${escapeHtml(item.source_file)}
      </button>
    </li>
  `).join("");
  elFileTabContent.innerHTML = results.map((item, index) => `
    <div class="tab-pane fade ${index === 0 ? "show active" : ""}" id="file-pane-${index}" role="tabpanel">
      ${renderFilePanel(item, index)}
    </div>
  `).join("");
  results.forEach((item, index) => renderFilePanelPlots(item, index));
}

function renderFilePanel(item, fileIndex) {
  const warningLines = item.warnings.map((warning) => formatWarningLine(warning)).join("");
  const outputLinks = Object.entries(item.outputs).map(([filename]) => `
    <a href="#" class="download-link" data-file-index="${fileIndex}" data-output-name="${escapeHtml(filename)}">
      <span>${escapeHtml(filename)}</span><span class="download-meta">CSV</span>
    </a>
  `).join("");
  return `
    <div class="file-panel">
      <div class="file-summary-grid">
        <div class="file-summary-block"><div class="file-summary-title">Input file</div><div class="file-summary-main">${escapeHtml(item.source_file)}</div></div>
        <div class="file-summary-block"><div class="file-summary-title">Rows with outliers</div><div class="file-summary-main">${item.summary.with_outliers_rows}</div></div>
        <div class="file-summary-block"><div class="file-summary-title">Rows without outliers</div><div class="file-summary-main">${item.summary.without_outliers_rows}</div></div>
        <div class="file-summary-block"><div class="file-summary-title">Removed outliers</div><div class="file-summary-main">${item.summary.removed_outliers_rows}</div></div>
        <div class="file-summary-block"><div class="file-summary-title">Warnings</div><div class="file-summary-main">${item.summary.warning_count}</div></div>
      </div>
      <div class="download-box"><div class="file-summary-title">Download per-file outputs</div><div class="download-files">${outputLinks}</div></div>
      <div class="warnings-box"><div class="file-summary-title">Warnings and processing notes</div><div class="warnings-list">${warningLines || '<div class="warning-line">No warnings.</div>'}</div></div>
      <div class="md-card" style="padding:18px;">
        <ul class="nav nav-tabs branch-tabs" role="tablist">
          <li class="nav-item"><button class="nav-link active" data-bs-toggle="tab" data-bs-target="#branch-with-${fileIndex}" type="button">with_outliers</button></li>
          <li class="nav-item"><button class="nav-link" data-bs-toggle="tab" data-bs-target="#branch-without-${fileIndex}" type="button">without_outliers</button></li>
        </ul>
        <div class="tab-content" style="padding-top:20px;">
          <div class="tab-pane fade show active" id="branch-with-${fileIndex}">${branchPlotsMarkup(fileIndex, "with_outliers")}</div>
          <div class="tab-pane fade" id="branch-without-${fileIndex}">${branchPlotsMarkup(fileIndex, "without_outliers")}</div>
        </div>
      </div>
    </div>
  `;
}

function branchPlotsMarkup(fileIndex, branch) {
  const baseId = `${branch}-${fileIndex}`;
  return `
    <div class="plots-grid">
      <div class="plot-card"><div class="plot-card-header"><div><h3 class="plot-title">Raw distribution by condition</h3><p class="plot-subtitle">Each point is one embryo. DMSO is always shown first.</p></div></div><div class="plot-shell"><div id="plot-raw-${baseId}" class="plotly-host"></div><div id="legend-raw-${baseId}" class="legend-stack"></div></div></div>
      <div class="plot-card"><div class="plot-card-header"><div><h3 class="plot-title">Raw median intensity across dates</h3><p class="plot-subtitle">Daily medians by condition before normalization.</p></div></div><div id="plot-rawline-${baseId}" class="plotly-host"></div></div>
      <div class="plot-card"><div class="plot-card-header"><div><h3 class="plot-title">DMSO anchor shift across dates</h3><p class="plot-subtitle">Internal reference used to correct day-to-day acquisition shifts.</p></div></div><div class="plot-shell"><div id="plot-dmso-${baseId}" class="plotly-host"></div><div id="legend-dmso-${baseId}" class="legend-stack"></div></div></div>
      <div class="plot-card"><div class="plot-card-header"><div><h3 class="plot-title">Normalized distribution by condition</h3><p class="plot-subtitle">Log2 ratios relative to the DMSO median from the same date.</p></div></div><div class="plot-shell"><div id="plot-norm-${baseId}" class="plotly-host"></div><div id="legend-norm-${baseId}" class="legend-stack"></div></div></div>
      <div class="plot-card"><div class="plot-card-header"><div><h3 class="plot-title">Variation reduction</h3><p class="plot-subtitle">Across-date CV of daily medians before and after normalization.</p></div></div><div id="plot-var-${baseId}" class="plotly-host"></div></div>
    </div>
  `;
}

function renderFilePanelPlots(item, fileIndex) {
  ["with_outliers", "without_outliers"].forEach((branch) => {
    const branchData = item.branches[branch];
    renderCategoricalScatterChart(`plot-raw-${branch}-${fileIndex}`, `legend-raw-${branch}-${fileIndex}`, branchData.long_rows, "intensity", "Raw intensity", false);
    renderRawMedianLineChart(`plot-rawline-${branch}-${fileIndex}`, branchData.summary_rows);
    renderDmsoAnchorChart(`plot-dmso-${branch}-${fileIndex}`, `legend-dmso-${branch}-${fileIndex}`, branchData.long_rows, branchData.dmso_rows);
    renderCategoricalScatterChart(`plot-norm-${branch}-${fileIndex}`, `legend-norm-${branch}-${fileIndex}`, branchData.long_rows.filter((row) => row.log2fc_vs_dmso !== null), "log2fc_vs_dmso", "log2(intensity / DMSO median)", true);
    renderVariationChart(`plot-var-${branch}-${fileIndex}`, branchData.variation_rows);
  });
  document.querySelectorAll(`[data-file-index="${fileIndex}"][data-output-name]`).forEach((link) => {
    link.addEventListener("click", (e) => {
      e.preventDefault();
      const outputName = link.dataset.outputName;
      downloadTextFile(outputName, item.outputs[outputName]);
    });
  });
}

function conditionOrder(rows) {
  const entries = [];
  const seen = new Set();
  rows.forEach((row) => {
    const key = row.condition_key;
    if (seen.has(key)) return;
    seen.add(key);
    entries.push({ key, label: row.condition_original });
  });
  entries.sort((a, b) => a.key.localeCompare(b.key));
  const dmso = entries.filter((item) => item.label.trim().toLowerCase() === "dmso");
  const others = entries.filter((item) => item.label.trim().toLowerCase() !== "dmso");
  return [...dmso, ...others].map((item) => item.label);
}

function buildColorMap(labels) {
  const map = {};
  labels.forEach((label, index) => { map[label] = SET2[index % SET2.length]; });
  return map;
}

function buildDateMarkerMap(dateLabels) {
  const map = {};
  dateLabels.forEach((label, index) => { map[label] = MARKERS[index % MARKERS.length]; });
  return map;
}

function buildPlotLegendHtml(colorMap, markerMap) {
  const conditionHtml = `<div class="legend-box"><div class="legend-title">Condition</div><div class="legend-items">${Object.entries(colorMap).map(([label, color]) => `<div class="legend-item"><span class="legend-swatch" style="background:${color};"></span><span>${escapeHtml(label)}</span></div>`).join("")}</div></div>`;
  const dateHtml = `<div class="legend-box"><div class="legend-title">Date</div><div class="legend-items">${Object.entries(markerMap).map(([label, marker]) => `<div class="legend-item"><span class="legend-symbol">${markerSymbolLabel(marker)}</span><span>${escapeHtml(label)}</span></div>`).join("")}</div></div>`;
  return conditionHtml + dateHtml;
}

function markerSymbolLabel(marker) {
  return { circle: "●", square: "■", "triangle-up": "▲", diamond: "◆", cross: "✚", x: "✖", "triangle-down": "▼", "triangle-left": "◀", "triangle-right": "▶", star: "★" }[marker] || "•";
}

function renderCategoricalScatterChart(plotId, legendId, rows, yKey, yTitle, drawZeroLine) {
  const plotEl = document.getElementById(plotId);
  const legendEl = document.getElementById(legendId);
  if (!plotEl || !legendEl) return;
  if (!rows.length) {
    plotEl.innerHTML = "<p>No data available for this branch.</p>";
    legendEl.innerHTML = "";
    return;
  }
  const order = conditionOrder(rows);
  const xPositions = Object.fromEntries(order.map((label, index) => [label, index]));
  const colorMap = buildColorMap(order);
  const dateLabels = [...new Set(rows.map((row) => row.acquisition_date))].sort();
  const markerMap = buildDateMarkerMap(dateLabels);
  legendEl.innerHTML = buildPlotLegendHtml(colorMap, markerMap);
  const traces = [];
  dateLabels.forEach((dateLabel) => {
    order.forEach((conditionLabel) => {
      const subset = rows.filter((row) => row.acquisition_date === dateLabel && row.condition_original === conditionLabel && row[yKey] !== null);
      if (!subset.length) return;
      const x = subset.map((_, idx) => xPositions[conditionLabel] + (((idx % 5) - 2) * 0.05));
      traces.push({
        type: "scatter", mode: "markers", x, y: subset.map((row) => row[yKey]),
        text: subset.map((row) => `${row.condition_original}<br>${row.acquisition_date}<br>${yTitle}: ${formatNumber(row[yKey])}`),
        hovertemplate: "%{text}<extra></extra>",
        marker: { size: 11, color: colorMap[conditionLabel], symbol: markerMap[dateLabel], line: { color: "#2d2d2d", width: 1.3 }, opacity: 0.9 },
        showlegend: false,
      });
    });
  });
  const shapes = [];
  order.forEach((conditionLabel) => {
    const values = rows.filter((row) => row.condition_original === conditionLabel && row[yKey] !== null).map((row) => row[yKey]);
    if (!values.length) return;
    const med = median(values);
    const xpos = xPositions[conditionLabel];
    shapes.push({ type: "line", x0: xpos - 0.23, x1: xpos + 0.23, y0: med, y1: med, line: { color: "#000000", width: 3 } });
  });
  if (drawZeroLine) shapes.push({ type: "line", x0: -0.5, x1: order.length - 0.5, y0: 0, y1: 0, line: { color: "#4c97d8", width: 2, dash: "dash" } });
  Plotly.newPlot(plotEl, traces, {
    margin: { l: 70, r: 10, t: 10, b: 110 }, paper_bgcolor: "white", plot_bgcolor: "white", showlegend: false,
    xaxis: { tickmode: "array", tickvals: order.map((_, idx) => idx), ticktext: order.map(wrapLabel), tickangle: -30, range: [-0.6, order.length - 0.4], title: "Condition", zeroline: false },
    yaxis: { title: yTitle, gridcolor: "#e7ebf0", zeroline: false }, shapes,
  }, { responsive: true, displaylogo: false });
}

function renderRawMedianLineChart(plotId, summaryRows) {
  const plotEl = document.getElementById(plotId);
  if (!plotEl) return;
  const order = conditionOrder(summaryRows);
  const colorMap = buildColorMap(order);
  const traces = order.map((conditionLabel) => {
    const subset = summaryRows.filter((row) => row.condition_original === conditionLabel).sort((a, b) => a.acquisition_date.localeCompare(b.acquisition_date));
    return { type: "scatter", mode: "lines+markers", x: subset.map((row) => row.acquisition_date), y: subset.map((row) => row.raw_median), name: conditionLabel, line: { color: colorMap[conditionLabel], width: 3 }, marker: { color: colorMap[conditionLabel], size: 8 } };
  });
  Plotly.newPlot(plotEl, traces, {
    margin: { l: 70, r: 15, t: 10, b: 70 }, paper_bgcolor: "white", plot_bgcolor: "white",
    legend: { orientation: "h", y: -0.25 }, xaxis: { title: "Acquisition date" }, yaxis: { title: "Median raw intensity", gridcolor: "#e7ebf0" },
  }, { responsive: true, displaylogo: false });
}

function renderDmsoAnchorChart(plotId, legendId, longRows, dmsoRows) {
  const plotEl = document.getElementById(plotId);
  const legendEl = document.getElementById(legendId);
  if (!plotEl || !legendEl) return;
  const dmsoPoints = longRows.filter((row) => row.is_dmso_condition);
  const order = dmsoRows.map((row) => row.acquisition_date);
  const xPositions = Object.fromEntries(order.map((label, index) => [label, index]));
  const colorMap = buildColorMap(order);
  const markerMap = buildDateMarkerMap(order);
  legendEl.innerHTML = `<div class="legend-box"><div class="legend-title">Date</div><div class="legend-items">${order.map((dateLabel) => `<div class="legend-item"><span class="legend-swatch" style="background:${colorMap[dateLabel]};"></span><span>${escapeHtml(dateLabel)}</span></div>`).join("")}</div></div>`;
  const traces = [];
  order.forEach((dateLabel) => {
    const subset = dmsoPoints.filter((row) => row.acquisition_date === dateLabel);
    const x = subset.map((_, idx) => xPositions[dateLabel] + (((idx % 5) - 2) * 0.04));
    traces.push({ type: "scatter", mode: "markers", x, y: subset.map((row) => row.intensity), marker: { size: 11, color: colorMap[dateLabel], symbol: markerMap[dateLabel], line: { color: "#2d2d2d", width: 1.3 } }, showlegend: false });
  });
  traces.push({ type: "scatter", mode: "lines+markers", x: order.map((dateLabel) => xPositions[dateLabel]), y: dmsoRows.map((row) => row.dmso_median), line: { color: "#000000", width: 2 }, marker: { color: "#000000", size: 8 }, showlegend: false });
  Plotly.newPlot(plotEl, traces, {
    margin: { l: 70, r: 10, t: 10, b: 90 }, paper_bgcolor: "white", plot_bgcolor: "white", showlegend: false,
    xaxis: { tickmode: "array", tickvals: order.map((_, idx) => idx), ticktext: order, tickangle: -25, title: "Acquisition date", range: [-0.5, order.length - 0.5] },
    yaxis: { title: "Raw DMSO intensity", gridcolor: "#e7ebf0" },
  }, { responsive: true, displaylogo: false });
}

function renderVariationChart(plotId, variationRows) {
  const plotEl = document.getElementById(plotId);
  if (!plotEl) return;
  const usable = variationRows.filter((row) => row.raw_daily_median_cv !== null && row.normalized_daily_median_cv !== null);
  if (!usable.length) {
    plotEl.innerHTML = "<p>No sufficient across-date data for this summary.</p>";
    return;
  }
  const order = usable.slice().sort((a, b) => (b.raw_daily_median_cv || 0) - (a.raw_daily_median_cv || 0)).map((row) => row.condition_original);
  const colorMap = buildColorMap(order);
  const traces = usable.flatMap((row) => {
    const y = row.condition_original;
    return [
      { type: "scatter", mode: "lines", x: [row.raw_daily_median_cv, row.normalized_daily_median_cv], y: [y, y], line: { color: "#b4bcc7", width: 2 }, hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "markers", x: [row.raw_daily_median_cv], y: [y], marker: { color: colorMap[y], size: 11, symbol: "circle", line: { color: "#2d2d2d", width: 1 } }, showlegend: false },
      { type: "scatter", mode: "markers", x: [row.normalized_daily_median_cv], y: [y], marker: { color: colorMap[y], size: 11, symbol: "square", line: { color: "#2d2d2d", width: 1 } }, showlegend: false },
    ];
  });
  Plotly.newPlot(plotEl, traces, {
    margin: { l: 120, r: 30, t: 10, b: 70 }, paper_bgcolor: "white", plot_bgcolor: "white",
    xaxis: { title: "Coefficient of variation of daily medians", gridcolor: "#e7ebf0" },
    yaxis: { categoryorder: "array", categoryarray: order },
    annotations: [{ xref: "paper", yref: "paper", x: 0.02, y: 1.08, text: "● Raw   ■ Normalized", showarrow: false, font: { size: 12, color: "#4d4d4d" } }],
  }, { responsive: true, displaylogo: false });
}

async function downloadAllOutputsZip() {
  if (!analysisResults) return;
  const zip = new JSZip();
  for (const result of analysisResults.results) {
    const folder = zip.folder(result.source_file.replace(/\.csv$/i, ""));
    Object.entries(result.outputs).forEach(([filename, content]) => {
      folder.file(filename, content || "");
    });
  }
  zip.file("run_summary.json", JSON.stringify(analysisResults.run_summary, null, 2));
  zip.file("ignored_files.json", JSON.stringify(analysisResults.ignored_files || [], null, 2));
  const blob = await zip.generateAsync({ type: "blob" });
  saveAs(blob, "zebrafish-ros-normalizer-online-results.zip");
}

function downloadTextFile(filename, content) {
  const blob = new Blob([content], { type: "text/csv;charset=utf-8" });
  saveAs(blob, filename);
}

function formatWarningLine(warning) {
  const pieces = [warning.level?.toUpperCase(), warning.code, warning.message];
  if (warning.source_file) pieces.push(`file=${warning.source_file}`);
  if (warning.row_number) pieces.push(`row=${warning.row_number}`);
  if (warning.column_name) pieces.push(`column=${warning.column_name}`);
  if (warning.date_raw) pieces.push(`date=${warning.date_raw}`);
  if (warning.value) pieces.push(`value=${warning.value}`);
  return `<div class="warning-line">${escapeHtml(pieces.filter(Boolean).join(" | "))}</div>`;
}

function deepConvert(value) {
  if (Array.isArray(value)) return value.map(deepConvert);
  if (value && typeof value === "object") {
    if (typeof value.toJs === "function") return deepConvert(value.toJs({ dict_converter: Object.fromEntries }));
    const out = {};
    for (const [key, val] of Object.entries(value)) out[key] = deepConvert(val);
    return out;
  }
  return value;
}

function median(values) {
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

function wrapLabel(label) {
  const text = String(label);
  if (text.length <= 10) return text;
  const parts = text.split(/\s+/);
  if (parts.length === 1) return text;
  const lines = [];
  let current = "";
  parts.forEach((part) => {
    if ((`${current} ${part}`).trim().length > 10 && current) {
      lines.push(current);
      current = part;
    } else {
      current = `${current} ${part}`.trim();
    }
  });
  if (current) lines.push(current);
  return lines.join("<br>");
}

function formatBytes(bytes) {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
}

function formatNumber(value) {
  if (value === null || value === undefined || Number.isNaN(value)) return "NA";
  return Number(value).toFixed(3);
}

function escapeHtml(text) {
  return String(text).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;").replace(/'/g, "&#039;");
}

function show(el) { el.classList.remove("d-none"); }
function hide(el) { el.classList.add("d-none"); }
function sleep(ms) { return new Promise((resolve) => setTimeout(resolve, ms)); }
