/*! © 2024 Nourhan Abdelfattah — scSidekick R package
 *  Bundled for use via scSidekick::ChartBuilder() only.
 *  Not licensed for standalone redistribution.
 */

// ═══════════════════════════════════════════════════════════════════
//  FILE UPLOAD  (depends on: config.js, ui_controls.js)
//  NOTE: upload UI is not present in the scSidekick bundle —
//  data always arrives via NOURKIT_PAYLOAD (injected by ChartBuilder()).
// ═══════════════════════════════════════════════════════════════════
const dropZone  = document.getElementById("dropZone");
const fileInput = document.getElementById("fileInput");

dropZone.addEventListener("dragover",  e => { e.preventDefault(); dropZone.classList.add("drag-over"); });
dropZone.addEventListener("dragleave", () => dropZone.classList.remove("drag-over"));
dropZone.addEventListener("drop",      e => { e.preventDefault(); dropZone.classList.remove("drag-over"); handleFile(e.dataTransfer.files[0]); });
fileInput.addEventListener("change",   e => handleFile(e.target.files[0]));

function handleFile(file) {
  if (!file) return;
  const ext = file.name.split(".").pop().toLowerCase();
  showSpin("Reading file…");
  if (ext === "csv") parseCSV(file);
  else if (["xlsx","xls"].includes(ext)) parseExcel(file);
  else { hideSpin(); alert("Please upload a .csv, .xlsx, or .xls file."); }
}

function parseCSV(file) {
  Papa.parse(file, {
    header: true, skipEmptyLines: true, dynamicTyping: false,
    complete: r => { loadData(r.data, file.name); hideSpin(); },
    error:    e => { hideSpin(); alert("CSV error: " + e.message); },
  });
}

function parseExcel(file) {
  const reader = new FileReader();
  reader.onload = e => {
    try {
      const wb   = XLSX.read(e.target.result, { type:"array" });
      const ws   = wb.Sheets[wb.SheetNames[0]];
      const rows = XLSX.utils.sheet_to_json(ws, { defval:"" });
      loadData(rows, file.name);
    } catch(err) { alert("Excel error: " + err.message); }
    hideSpin();
  };
  reader.readAsArrayBuffer(file);
}

// ═══════════════════════════════════════════════════════════════════
//  DATA LOADING + TYPE DETECTION
// ═══════════════════════════════════════════════════════════════════
function loadData(rows, filename) {
  if (!rows.length) { alert("File appears to be empty."); return; }

  // ── Fix empty / __EMPTY column headers ──────────────────────────
  const rawKeys = Object.keys(rows[0]);
  const keyMap  = {};
  let emptyIdx  = 0;
  rawKeys.forEach(k => {
    const trimmed = String(k).trim();
    if (!trimmed || trimmed === "__EMPTY" || /^__EMPTY_?\d*$/.test(trimmed)) {
      emptyIdx++;
      keyMap[k] = emptyIdx === 1 ? "Feature" : `Feature_${emptyIdx}`;
    } else {
      keyMap[k] = k;
    }
  });
  const needsRemap = rawKeys.some(k => keyMap[k] !== k);
  if (needsRemap) {
    rows = rows.map(row => {
      const nr = {};
      rawKeys.forEach(k => { nr[keyMap[k]] = row[k]; });
      return nr;
    });
  }

  // ── Sample / save ────────────────────────────────────────────────
  RAW_DATA = rows;
  DATA = rows.length > MAX_ROWS
    ? [...rows].sort(() => Math.random() - .5).slice(0, MAX_ROWS)
    : rows;

  // Reset explorer state on new load
  COL_OVERRIDES   = {};
  COL_RENAMES     = {};
  ROW_FILTER_TEXT = "";
  _explorerPage   = 0;

  const keys = Object.keys(DATA[0]);
  COLS = keys.map(name => {
    const vals     = DATA.map(r => r[name]).filter(v => v !== null && v !== "" && v !== undefined);
    const type     = COL_OVERRIDES[name] || detectType(vals);
    const nums     = type === "num" ? vals.map(v => parseFloat(v)).filter(v => !isNaN(v)) : [];
    const uniqVals = [...new Set(vals)];
    return { name, type, uniq: uniqVals.length,
             min: nums.length ? Math.min(...nums) : null,
             max: nums.length ? Math.max(...nums) : null };
  }).filter(c => c.type !== "empty");

  _pendingFilename = filename;

  // Always route through the preview/setup screen so the user can inspect
  // the data, fix column types, rename columns, and confirm the format
  // before plotting — regardless of whether it looks long or wide.
  FORMAT = detectWideFormat() ? "wide" : "long";
  showSetupScreen();
}

function detectType(vals) {
  if (!vals.length) return "empty";
  const n    = vals.length;
  const numN = vals.filter(v => v !== "" && !isNaN(parseFloat(v)) && isFinite(v)).length;
  if (numN / n >= 0.85) return "num";
  const dateN = vals.filter(v => { try { return !!v && isNaN(+v) && !isNaN(Date.parse(v)); } catch(e){return false;} }).length;
  if (dateN / n >= 0.8) return "date";
  return "cat";
}

// ── Structure-based wide format detection ────────────────────────────
//
// Wide format (feature × sample matrix): one row per feature/entity,
//   numeric columns represent samples, at least one categorical column
//   is a unique row identifier.
//
// Long format (tidy): rows repeat across conditions, one "value" column
//   holds the measurement, categorical columns are grouping variables.
//
// The toggle on the setup screen lets the user override the guess.

function detectWideFormat() {
  if (!COLS.length || !DATA.length) return false;
  const n        = DATA.length;
  const numCols  = COLS.filter(c => c.type === "num");
  const catCols  = COLS.filter(c => c.type !== "num");
  const numN     = numCols.length;

  if (numN < 2) return false;

  // ── Long-format signal 1: a "value" column with a measurement name ──
  const VALUE_RE = /^(values?|n|counts?|means?|avg|averages?|expression|levels?|concentrations?|measurements?|scores?|signals?|results?|intensit\w*|abundance|amounts?|quantit\w*|log2?fc?|fold.?change)$/i;
  const hasValueName = numCols.some(c => VALUE_RE.test(c.name.trim()));

  // ── Long-format signal 2: rows repeat across grouping variables ──────
  // A "grouping variable" is a categorical column with low cardinality
  // (< 30 % unique values) — the hallmark of a tidy grouping factor.
  const groupingCols = catCols.filter(c => c.uniq / n < 0.30);
  // Only flag repetition when there is NO high-uniqueness row-ID column
  // (wide format also has grouping-like annotation columns).
  const hasRowId = catCols.some(c => c.uniq / n >= 0.85);
  let rowsRepeat = false;
  if (groupingCols.length > 0 && !hasRowId) {
    const gcNames      = groupingCols.map(c => c.name);
    const uniqueCombos = new Set(
      DATA.map(r => gcNames.map(cn => String(r[cn] ?? "")).join("\x00"))
    );
    rowsRepeat = uniqueCombos.size / n < 0.50;
  }

  // Clear long-format evidence → not wide
  if (hasValueName && rowsRepeat)  return false;
  if (hasValueName && numN <= 3)   return false;
  if (rowsRepeat   && numN <= 4)   return false;

  // ── Wide-format fast path: many numeric columns ──────────────────────
  if (numN >= 10 && !hasValueName) return true;

  // ── Wide-format structural evidence ─────────────────────────────────
  // Requires BOTH a unique row identifier AND a sample naming pattern.
  const hasSampleNames = _hasSampleNamingPattern(numCols.map(c => c.name));
  return hasRowId && hasSampleNames;
}

// Returns true when numeric column names share a systematic group+replicate
// pattern, e.g. "WT.1","WT.2","KO.1","KO.2" or "Sample1"…"SampleN".
function _hasSampleNamingPattern(names) {
  if (names.length < 3) return false;

  // Delimiter-based grouping: "WT.1", "KO_2", "Ctrl-A"
  for (const d of [".", "_", "-", " "]) {
    const withDelim = names.filter(n => n.includes(d));
    if (withDelim.length < names.length * 0.5) continue;
    const prefixes  = withDelim.map(n => n.split(d)[0]);
    const uniqPre   = new Set(prefixes);
    // Multiple groups (≥2) with fewer groups than columns  OR  single prefix ≥4 cols
    if (uniqPre.size >= 2 && uniqPre.size <= names.length * 0.65) return true;
    if (uniqPre.size === 1 && names.length >= 4)                  return true;
  }

  // Numeric-suffix grouping: "WT1","KO2" or "Mouse1"…"Mouse8"
  const numSuffix = names.filter(n => /\d+$/.test(n));
  if (numSuffix.length >= names.length * 0.70) {
    const prefixes = numSuffix.map(n => n.replace(/\d+$/, ""));
    const uniqPre  = new Set(prefixes);
    if (uniqPre.size >= 2 && uniqPre.size <= names.length * 0.65) return true;
    if (uniqPre.size === 1 && names.length >= 4)                  return true;
  }

  return false;
}

function suggestCharts(cols) {
  const nums  = cols.filter(c => c.type === "num");
  const cats  = cols.filter(c => c.type === "cat");
  const dates = cols.filter(c => c.type === "date");
  const out   = [];
  if (cats.length >= 1 && nums.length >= 1) out.push("bar", "box", "violin");
  if (cats.length >= 2 && nums.length >= 1) out.push("heatmap");
  if (nums.length >= 2)                     out.push("scatter");
  if (dates.length >= 1 && nums.length >= 1) out.push("line");
  if (nums.length >= 1 && cats.length === 0) out.push("histogram");
  if (cats.length >= 2)                     out.push("alluvial");
  if (!out.length) out.push("bar");
  return [...new Set(out)];
}

// ═══════════════════════════════════════════════════════════════════
//  UNIVERSAL SETUP / PREVIEW SCREEN
//  Shown for every file load. Lets the user fix column types,
//  rename columns, choose Long vs Wide, configure Wide settings,
//  then Visualize.
// ═══════════════════════════════════════════════════════════════════
function showSetupScreen() {
  // Sync format toggle buttons inside the setup screen
  document.querySelectorAll(".setup-fmt-btn").forEach(b =>
    b.classList.toggle("on", b.dataset.fmt === FORMAT)
  );

  // Build / refresh the editable data-preview table
  buildSetupPreview();

  // Show / hide long vs wide sections and update hints
  _applySetupFormat(FORMAT);

  // Transition screens
  document.getElementById("uploadScreen").style.display   = "none";
  document.getElementById("appScreen").classList.remove("visible");
  document.getElementById("setupScreen").classList.add("visible");
}

// Called by the Long / Wide toggle buttons inside the setup screen
function switchSetupFormat(fmt) {
  if (fmt === FORMAT) return;
  FORMAT = fmt;
  document.querySelectorAll(".setup-fmt-btn").forEach(b =>
    b.classList.toggle("on", b.dataset.fmt === fmt)
  );
  _applySetupFormat(fmt);
}

// Show/hide sections + update text for the current format choice
function _applySetupFormat(fmt) {
  const isWide = fmt === "wide";

  document.getElementById("setupWideSection").style.display = isWide ? "" : "none";
  document.getElementById("setupLongSection").style.display = isWide ? "none" : "";

  const confirmBtn = document.getElementById("setupConfirmBtn");
  if (confirmBtn) confirmBtn.textContent = isWide ? "Confirm & Visualize →" : "Visualize →";

  const hint = document.getElementById("setupFormatHint");
  if (hint) hint.textContent = isWide
    ? "Wide / matrix: one row per feature, columns are samples. Best for heatmaps."
    : "Long / tidy: one row per observation. Best for bar, box, scatter, and line charts.";

  if (isWide) _initWideSetupUI();
}

// Initialise the wide-format meta-column checkboxes + delimiter UI
function _initWideSetupUI() {
  const allCols = Object.keys(DATA[0]);
  const checkDiv = document.getElementById("metaCheckboxes");
  checkDiv.innerHTML = "";
  allCols.forEach(col => {
    const colInfo = COLS.find(c => c.name === col);
    const isMeta  = col === allCols[0] || !colInfo || colInfo.type !== "num";
    const div = document.createElement("div");
    div.className = "meta-check" + (isMeta ? " on" : "");
    div.dataset.col = col;
    div.onclick = () => { div.classList.toggle("on"); refreshSetupDelim(); };
    div.innerHTML = `<span>${esc(col)}</span>`;
    checkDiv.appendChild(div);
  });
  const numColNames = COLS.filter(c => c.type === "num").map(c => c.name);
  const bestDelim   = autoDetectDelimiter(numColNames);
  WS.delimiter = bestDelim;
  document.querySelectorAll(".delim-pill").forEach(p =>
    p.classList.toggle("on", p.dataset.delim === bestDelim)
  );
  refreshSetupDelim();
}

// ── Editable preview table ────────────────────────────────────────────
function buildSetupPreview() {
  const container = document.getElementById("setupPreviewTable");
  if (!container || !DATA.length || !COLS.length) return;

  const typeCycle = { cat:"num", num:"date", date:"cat" };
  const typeCls   = { num:"tb-num", cat:"tb-cat", date:"tb-date" };
  const showCols  = COLS.slice(0, 10);
  const showRows  = DATA.slice(0, 8);

  const thCells = showCols.map(col => {
    const tb = typeCls[col.type] || "tb-any";
    return `<th class="sp-th">
      <div class="sp-col-head">
        <span class="sp-col-name" title="Double-click to rename"
              ondblclick="promptSetupRenameCol('${esc(col.name)}')">${esc(col.name)}</span>
        <span class="type-badge ${tb}" style="cursor:pointer" title="Click to change type"
              onclick="setupOverrideColType('${esc(col.name)}','${typeCycle[col.type]}')">${col.type}</span>
      </div>
    </th>`;
  }).join("");

  const trRows = showRows.map((row, ri) => {
    const tds = showCols.map(col => {
      const v       = row[col.name];
      const isEmpty = v === null || v === "" || v === undefined ||
                      (col.type === "num" && isNaN(parseFloat(v)));
      const display = isEmpty ? "" : String(v).slice(0, 30);
      return `<td class="${isEmpty ? "sp-null" : ""}">${esc(display)}</td>`;
    }).join("");
    return `<tr><td class="sp-rn">${ri + 1}</td>${tds}</tr>`;
  }).join("");

  const moreColsRow = COLS.length > 10
    ? `<tr><td class="sp-rn"></td>
         <td colspan="${showCols.length}" class="sp-more">
           … and ${COLS.length - 10} more column${COLS.length - 10 > 1 ? "s" : ""}
         </td></tr>`
    : "";

  const footerText = DATA.length > 8
    ? `${DATA.length.toLocaleString()} rows · ${COLS.length} columns — showing first 8 rows`
    : `${DATA.length.toLocaleString()} row${DATA.length > 1 ? "s" : ""} · ${COLS.length} column${COLS.length > 1 ? "s" : ""}`;

  container.innerHTML = `
    <div class="sp-wrap">
      <table class="sp-table">
        <thead><tr><th class="sp-rn">#</th>${thCells}</tr></thead>
        <tbody>${trRows}${moreColsRow}</tbody>
      </table>
    </div>
    <p class="sp-footer">${footerText} · <em>Click a type badge to change it · Double-click a column name to rename it</em></p>`;
}

// Change a column's type in the preview (updates COLS, rebuilds table + wide UI)
function setupOverrideColType(colName, newType) {
  COL_OVERRIDES[colName] = newType;
  const idx = COLS.findIndex(c => c.name === colName);
  if (idx >= 0) {
    const vals = DATA.map(r => r[colName]).filter(v => v !== null && v !== "" && v !== undefined);
    const nums = newType === "num" ? vals.map(v => parseFloat(v)).filter(v => !isNaN(v)) : [];
    COLS[idx] = { ...COLS[idx], type: newType,
      min: nums.length ? Math.min(...nums) : null,
      max: nums.length ? Math.max(...nums) : null };
  }
  buildSetupPreview();
  if (FORMAT === "wide") _initWideSetupUI();
}

// Rename a column in DATA, RAW_DATA and COLS, then rebuild the preview
// BUG FIX A: surgically update the checkbox instead of rebuilding all checkboxes
function setupRenameCol(oldName, newName) {
  newName = newName.trim();
  if (!newName || newName === oldName) return;
  if (COLS.some(c => c.name === newName)) { alert(`Column "${newName}" already exists.`); return; }
  COL_RENAMES[oldName] = newName;
  DATA.forEach(row    => { row[newName] = row[oldName]; delete row[oldName]; });
  RAW_DATA.forEach(row => { row[newName] = row[oldName]; delete row[oldName]; });
  const idx = COLS.findIndex(c => c.name === oldName);
  if (idx >= 0) COLS[idx] = { ...COLS[idx], name: newName };
  buildSetupPreview();
  if (FORMAT === "wide") {
    // Surgically update the one renamed checkbox — preserves meta/sample assignments
    const checkDiv = document.getElementById("metaCheckboxes");
    const el = checkDiv?.querySelector(`[data-col="${CSS.escape(oldName)}"]`);
    if (el) {
      el.dataset.col = newName;
      const span = el.querySelector("span");
      if (span) span.textContent = newName;
    }
    buildSetupPreview();
    refreshSetupDelim();
  }
}

function promptSetupRenameCol(oldName) {
  const newName = prompt(`Rename "${oldName}" to:`, oldName);
  if (newName && newName.trim() !== oldName) setupRenameCol(oldName, newName.trim());
}

function autoDetectDelimiter(sampleCols) {
  const delims = [".", "_", "-", " "];
  let best = ".", bestScore = -1;
  for (const d of delims) {
    const groups  = sampleCols.map(c => c.split(d)[0]);
    const uniqueG = [...new Set(groups)];
    if (uniqueG.length < 2 || uniqueG.length >= sampleCols.length) continue;
    const minPerGroup = Math.min(...uniqueG.map(g => groups.filter(x => x === g).length));
    const score = uniqueG.length * minPerGroup;
    if (score > bestScore) { bestScore = score; best = d; }
  }
  return bestScore > 0 ? best : ".";
}

function setDelim(d) {
  WS.delimiter = d;
  document.querySelectorAll(".delim-pill").forEach(p =>
    p.classList.toggle("on", p.dataset.delim === d)
  );
  refreshSetupDelim();
}

function refreshSetupDelim() {
  const metaChecked = [...document.querySelectorAll(".meta-check.on")].map(el => el.dataset.col);
  const allCols     = Object.keys(DATA[0]);
  const sampleCols  = allCols.filter(c => !metaChecked.includes(c));

  const groups = parseGroupsFromDelim(sampleCols, WS.delimiter);

  const seen = new Set();
  const newOrder = [];
  sampleCols.forEach(sc => {
    const g = groups[sc];
    if (!seen.has(g)) { seen.add(g); newOrder.push(g); }
  });
  const preserved = WS.groupOrder.filter(g => newOrder.includes(g));
  const added     = newOrder.filter(g => !preserved.includes(g));
  WS.groupOrder   = [...preserved, ...added];

  buildDelimPreview(sampleCols, groups);

  // Hide group order panel when there's only one group (no grouping selected)
  const hasGroups = WS.groupOrder.length > 1 && WS.delimiter !== "";
  document.getElementById("groupOrderSection").style.display = hasGroups ? "" : "none";
  buildGroupOrderUI("setupGroupOrder");
}

// Problem 3.5: when delimiter is empty, all samples fall into one group called "All"
function parseGroupsFromDelim(sampleCols, delim) {
  const out = {};
  sampleCols.forEach(sc => {
    out[sc] = delim ? sc.split(delim)[0] : "All";
  });
  return out;
}

function buildDelimPreview(sampleCols, groups) {
  const el = document.getElementById("delimPreview");
  if (!sampleCols.length) { el.innerHTML = ""; return; }
  const show = sampleCols.slice(0, 8);
  el.innerHTML = `<div class="preview-wrap"><table class="preview-table">
    <tr><th>Column name</th><th>Detected group</th></tr>
    ${show.map(sc => `<tr><td>${esc(sc)}</td><td><b>${esc(groups[sc])}</b></td></tr>`).join("")}
    ${sampleCols.length > 8 ? `<tr><td colspan="2" style="color:#9aa5b4;font-style:italic">… and ${sampleCols.length - 8} more</td></tr>` : ""}
  </table></div>`;
}

function buildGroupOrderUI(containerId) {
  const container = document.getElementById(containerId);
  if (!container) return;
  const colors = getColors();
  container.innerHTML = "";
  WS.groupOrder.forEach((g, i) => {
    const item = document.createElement("div");
    item.className = "group-order-item";
    const color = colors[i % colors.length];
    item.innerHTML = `
      <div class="group-color-swatch" style="background:${color}"></div>
      <span class="group-item-name">${esc(g)}</span>
      <div class="move-btns">
        <button class="move-btn" onclick="moveGroup(-1,${i})" ${i === 0 ? "disabled" : ""}>▲</button>
        <button class="move-btn" onclick="moveGroup(1,${i})"  ${i === WS.groupOrder.length - 1 ? "disabled" : ""}>▼</button>
      </div>`;
    container.appendChild(item);
  });
}

function moveGroup(dir, idx) {
  const ni = idx + dir;
  if (ni < 0 || ni >= WS.groupOrder.length) return;
  [WS.groupOrder[idx], WS.groupOrder[ni]] = [WS.groupOrder[ni], WS.groupOrder[idx]];
  buildGroupOrderUI("setupGroupOrder");
  buildGroupOrderUI("sidebarGroupOrder");
  if (FORMAT === "wide" && ST.chartType === "heatmap") render();
}

function confirmSetup() {
  document.getElementById("setupScreen").classList.remove("visible");

  if (FORMAT === "long") {
    initApp(_pendingFilename);
    return;
  }

  // Wide format: finalize WS configuration
  WS.metaCols   = [...document.querySelectorAll(".meta-check.on")].map(el => el.dataset.col);
  const allCols = Object.keys(DATA[0]);
  WS.sampleCols = allCols.filter(c => !WS.metaCols.includes(c));

  if (!WS.sampleCols.length) {
    // Put the screen back so the user can fix it
    document.getElementById("setupScreen").classList.add("visible");
    alert("Please leave at least one column unchecked as sample/measurement data.");
    return;
  }
  // (no hard block on zero meta cols — renderWideHeatmap uses row index when rowLabelCol = null)

  WS.sampleGroups = parseGroupsFromDelim(WS.sampleCols, WS.delimiter);
  const uniqueGroups = [...new Set(Object.values(WS.sampleGroups))];
  ST.colSplit = uniqueGroups.length > 1;

  initApp(_pendingFilename);
}

// ═══════════════════════════════════════════════════════════════════
//  FORMAT SWITCHER  (Problem 3.2)
// ═══════════════════════════════════════════════════════════════════
function switchFormat(newFormat) {
  if (newFormat === FORMAT) return;
  if (newFormat === "wide") {
    FORMAT = "wide";
    // Reset wide state so setup screen feels fresh
    WS.metaCols = []; WS.sampleCols = []; WS.sampleGroups = {};
    WS.groupOrder = []; WS.delimiter = ".";
    WS.sampleMeta = {}; WS.metaFields = [];
    WS.rowFilter = null; WS.sampleFilter = null;
    document.getElementById("appScreen").classList.remove("visible");
    showSetupScreen();
  } else {
    FORMAT = "long";
    ST.colFilters = {};
    // Re-compute COLS from all columns and init as long format
    const keys = Object.keys(DATA[0]);
    COLS = keys.map(name => {
      const vals     = DATA.map(r => r[name]).filter(v => v !== null && v !== "" && v !== undefined);
      const type     = COL_OVERRIDES[name] || detectType(vals);
      const nums     = type === "num" ? vals.map(v => parseFloat(v)).filter(v => !isNaN(v)) : [];
      const uniqVals = [...new Set(vals)];
      return { name, type, uniq: uniqVals.length,
               min: nums.length ? Math.min(...nums) : null,
               max: nums.length ? Math.max(...nums) : null };
    }).filter(c => c.type !== "empty");
    initApp(_pendingFilename);
  }
  // Update format toggle button states
  document.querySelectorAll(".fmt-toggle-btn").forEach(b =>
    b.classList.toggle("on", b.dataset.fmt === FORMAT)
  );
}

// ═══════════════════════════════════════════════════════════════════
//  APP INIT
// ═══════════════════════════════════════════════════════════════════
function initApp(filename) {
  // Reset stats and legend whenever a new file is loaded
  WS.rowFilter = null; WS.sampleFilter = null;
  ST.colFilters = {};
  ST.groupColors = {};
  ST.heatScaleType = "seq";
  try {
    if (typeof clearStats === "function") clearStats();
  } catch(e) {}
  try {
    const figTa = document.getElementById("legendFigureTextarea");
    if (figTa) figTa.value = "";
    const statsTa = document.getElementById("legendStatsTextarea");
    if (statsTa) statsTa.value = "";
  } catch(e) {}

  if (FORMAT === "wide") {
    ST.suggested  = ["heatmap"];
    ST.chartType  = "heatmap";
    ST.log2       = false;
    ST.center     = false;
    document.getElementById("wideBadge").style.display = "";
  } else {
    ST.suggested  = suggestCharts(COLS);
    ST.chartType  = ST.suggested[0] || "bar";
    document.getElementById("wideBadge").style.display = "none";
  }

  document.getElementById("fileBadge").textContent = filename;
  document.getElementById("rowBadge").textContent  =
    FORMAT === "wide"
      ? `${DATA.length.toLocaleString()} features · ${WS.sampleCols.length} samples`
      : `${DATA.length.toLocaleString()} rows · ${COLS.length} cols`;

  // Update format toggle buttons
  document.querySelectorAll(".fmt-toggle-btn").forEach(b =>
    b.classList.toggle("on", b.dataset.fmt === FORMAT)
  );

  buildColPills();
  buildChartGrid();
  buildPalGrid();
  populateSelectors();

  if (FORMAT === "wide") {
    const sync = (id, val) => {
      const el = document.getElementById(id);
      if (el) el.classList.toggle("on", val);
    };
    sync("togLog2",     ST.log2);
    sync("togCenter",   ST.center);
    sync("togColSplit", ST.colSplit);
    buildGroupOrderUI("sidebarGroupOrder");
  } else {
    applyDefaults();
  }

  // Problem 2: purge old chart, then render after one browser paint
  Plotly.purge("plot");
  document.getElementById("appScreen").classList.add("visible");
  switchTab("chart");
  // Always show the legend panel immediately — content is populated after first render
  try {
    const lp = document.getElementById("legendPanel");
    if (lp) lp.style.display = "";
  } catch(e) {}
  requestAnimationFrame(() => render());
}

// ═══════════════════════════════════════════════════════════════════
//  FEATURE 1 — Back to Setup button
// ═══════════════════════════════════════════════════════════════════
function goBackToSetup() {
  document.getElementById("appScreen").classList.remove("visible");
  document.getElementById("setupScreen").classList.add("visible");
  buildSetupPreview();
  if (FORMAT === "wide") _initWideSetupUIPreservingState();
}

function _initWideSetupUIPreservingState() {
  const allCols = Object.keys(DATA[0]);
  const checkDiv = document.getElementById("metaCheckboxes");
  checkDiv.innerHTML = "";
  allCols.forEach(col => {
    const isMeta = WS.metaCols.includes(col);
    const div = document.createElement("div");
    div.className = "meta-check" + (isMeta ? " on" : "");
    div.dataset.col = col;
    div.onclick = () => { div.classList.toggle("on"); refreshSetupDelim(); };
    div.innerHTML = `<span>${esc(col)}</span>`;
    checkDiv.appendChild(div);
  });
  document.querySelectorAll(".delim-pill").forEach(p =>
    p.classList.toggle("on", p.dataset.delim === WS.delimiter)
  );
  // Show the metadata upload section again
  const sec = document.getElementById("setupWideSection");
  if (sec) sec.style.display = "";
  refreshSetupDelim();
  // Restore metadata upload status display if sampleMeta is loaded
  if (WS.metaFields?.length) {
    const matched = WS.sampleCols.filter(sc => WS.sampleMeta?.[sc]).length;
    const info = document.getElementById("metaUploadInfo");
    if (info) info.textContent =
      `✓ ${WS.metaFields.length} field${WS.metaFields.length !== 1 ? "s" : ""} · ${matched}/${WS.sampleCols.length} samples matched`;
    const status = document.getElementById("metaUploadStatus");
    if (status) status.style.display = "flex";
  }
}

// ═══════════════════════════════════════════════════════════════════
//  FEATURE 2 — Sample metadata CSV upload
// ═══════════════════════════════════════════════════════════════════
function loadSampleMetadata(event) {
  const file = event.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = e => {
    const parsed = Papa.parse(e.target.result.trim(), { header: true, skipEmptyLines: true });
    const rows = parsed.data;
    if (!rows.length) return;
    const headers = Object.keys(rows[0]);
    const sampleKey = headers[0];
    const fields    = headers.slice(1);

    WS.sampleMeta = {};
    rows.forEach(row => {
      const sn = String(row[sampleKey] ?? "");
      if (!sn) return;
      const meta = {};
      fields.forEach(f => { meta[f] = row[f]; });
      WS.sampleMeta[sn] = meta;
    });
    WS.metaFields = fields;

    // Derive the current sample columns from the live checkbox state
    // (WS.sampleCols is only finalized by confirmSetup; at upload time it may be empty)
    const currentMetaCols = [...document.querySelectorAll(".meta-check.on")].map(el => el.dataset.col);
    const allDataCols     = DATA.length ? Object.keys(DATA[0]) : [];
    const effectiveSampleCols = allDataCols.filter(c => !currentMetaCols.includes(c));

    // Override sampleGroups with the FIRST metadata field (primary group)
    if (fields.length) {
      const primary = fields[0];
      effectiveSampleCols.forEach(sc => {
        const m = WS.sampleMeta[sc];
        if (m && m[primary] != null) WS.sampleGroups[sc] = String(m[primary]);
      });
    }
    // Recompute groupOrder
    const seen = new Set(), newOrder = [];
    effectiveSampleCols.forEach(sc => {
      const g = WS.sampleGroups[sc] || sc;
      if (!seen.has(g)) { seen.add(g); newOrder.push(g); }
    });
    WS.groupOrder = newOrder;

    const matched = effectiveSampleCols.filter(sc => WS.sampleMeta[sc]).length;
    const info = document.getElementById("metaUploadInfo");
    if (info) info.textContent =
      `✓ ${fields.length} field${fields.length !== 1 ? "s" : ""} · ${matched}/${effectiveSampleCols.length} samples matched`;
    const status = document.getElementById("metaUploadStatus");
    if (status) status.style.display = "flex";
    const inp = document.getElementById("metaFileInput");
    if (inp) inp.value = "";

    buildGroupOrderUI("setupGroupOrder");
    if (FORMAT === "wide" && ST.chartType === "heatmap") render();
  };
  reader.readAsText(file);
}

function clearSampleMetadata() {
  WS.sampleMeta  = {};
  WS.metaFields  = [];
  // Re-derive groups from delimiter
  WS.sampleGroups = parseGroupsFromDelim(WS.sampleCols, WS.delimiter);
  const seen = new Set(), newOrder = [];
  WS.sampleCols.forEach(sc => {
    const g = WS.sampleGroups[sc] || sc;
    if (!seen.has(g)) { seen.add(g); newOrder.push(g); }
  });
  WS.groupOrder = newOrder;
  const status = document.getElementById("metaUploadStatus");
  if (status) status.style.display = "none";
  const inp = document.getElementById("metaFileInput");
  if (inp) inp.value = "";
  buildGroupOrderUI("setupGroupOrder");
  if (FORMAT === "wide" && ST.chartType === "heatmap") render();
}

// ═══════════════════════════════════════════════════════════════════
//  NOURKIT LAUNCHER — required entry point (scSidekick bundle)
// ═══════════════════════════════════════════════════════════════════
// This bundle REQUIRES a NOURKIT_PAYLOAD injected by scSidekick::ChartBuilder().
// If opened directly (no payload), the lock screen is shown and the app
// does not initialise.
window.addEventListener("load", () => {
  const lockScreen   = document.getElementById("lockScreen");
  const uploadScreen = document.getElementById("uploadScreen");

  if (typeof NOURKIT_PAYLOAD === "undefined" || !NOURKIT_PAYLOAD) {
    // No payload — show lock screen, suppress everything else.
    if (uploadScreen) uploadScreen.style.display = "none";
    if (lockScreen)   lockScreen.style.display   = "flex";
    return;
  }

  // Payload present — hide lock screen, proceed.
  if (lockScreen)   lockScreen.style.display   = "none";
  if (uploadScreen) uploadScreen.style.display = "none";

  const rows = NOURKIT_PAYLOAD.rows;
  if (!Array.isArray(rows) || !rows.length) {
    console.warn("NourKit payload present but contained no rows.");
    return;
  }
  try {
    loadData(rows, NOURKIT_PAYLOAD.filename || "nourkit_data.csv");
  } catch (e) {
    console.error("NourKit auto-load failed:", e);
  }
});
