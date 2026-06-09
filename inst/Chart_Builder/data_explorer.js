/*! © 2024–2026 Nourhan Abdelfattah — scSidekick R package
 *  Bundled for use via scSidekick::ChartBuilder() only.
 *  Not licensed for standalone redistribution.
 */
// ═══════════════════════════════════════════════════════════════════
//  DATA EXPLORER  (depends on: config.js, ui_controls.js)
// ═══════════════════════════════════════════════════════════════════

// ── Apply row filter and rebuild DATA from RAW_DATA ─────────────────
function applyRowFilter(text) {
  ROW_FILTER_TEXT = text.trim().toLowerCase();
  DATA = ROW_FILTER_TEXT
    ? RAW_DATA.filter(row =>
        Object.values(row).some(v => String(v).toLowerCase().includes(ROW_FILTER_TEXT))
      )
    : [...RAW_DATA];
  // Respect MAX_ROWS cap
  if (DATA.length > MAX_ROWS) DATA = DATA.slice(0, MAX_ROWS);
  _explorerPage = 0;
  renderDataTable();
  // Update row badge
  document.getElementById("rowBadge").textContent =
    FORMAT === "wide"
      ? `${DATA.length.toLocaleString()} features · ${WS.sampleCols.length} samples`
      : `${DATA.length.toLocaleString()} rows · ${COLS.length} cols`;
}

// ── Override a column's detected type ───────────────────────────────
function overrideColType(colName, newType) {
  COL_OVERRIDES[colName] = newType;
  // Recompute COLS entry
  const idx = COLS.findIndex(c => c.name === colName);
  if (idx >= 0) {
    const vals = DATA.map(r => r[colName]).filter(v => v !== null && v !== "" && v !== undefined);
    const nums = newType === "num" ? vals.map(v => parseFloat(v)).filter(v => !isNaN(v)) : [];
    COLS[idx] = { ...COLS[idx], type: newType,
      min: nums.length ? Math.min(...nums) : null,
      max: nums.length ? Math.max(...nums) : null };
  }
  buildColPills();
  buildChartGrid();
  buildPalGrid();
  populateSelectors();
  renderDataTable();
}

// ── Rename a column ──────────────────────────────────────────────────
function renameCol(oldName, newName) {
  newName = newName.trim();
  if (!newName || newName === oldName) return;
  // Check for conflicts
  if (COLS.some(c => c.name === newName)) {
    alert(`Column "${newName}" already exists.`); return;
  }
  COL_RENAMES[oldName] = newName;
  // Update DATA rows
  DATA.forEach(row => { row[newName] = row[oldName]; delete row[oldName]; });
  RAW_DATA.forEach(row => { row[newName] = row[oldName]; delete row[oldName]; });
  // Update COLS
  const idx = COLS.findIndex(c => c.name === oldName);
  if (idx >= 0) COLS[idx] = { ...COLS[idx], name: newName };
  // Update WS arrays if wide
  if (FORMAT === "wide") {
    const mi = WS.metaCols.indexOf(oldName);
    if (mi >= 0) WS.metaCols[mi] = newName;
    const si = WS.sampleCols.indexOf(oldName);
    if (si >= 0) {
      WS.sampleCols[si] = newName;
      if (WS.sampleGroups[oldName] !== undefined) {
        WS.sampleGroups[newName] = WS.sampleGroups[oldName];
        delete WS.sampleGroups[oldName];
      }
    }
  }
  buildColPills();
  populateSelectors();
  renderDataTable();
}

// ── Fill null / empty values in a column ────────────────────────────
function fillNulls(strategy) {
  const numCols = COLS.filter(c => c.type === "num").map(c => c.name);
  numCols.forEach(col => {
    const vals = DATA.map(r => parseFloat(r[col])).filter(v => !isNaN(v));
    if (!vals.length) return;
    let fillVal;
    if (strategy === "zero")   fillVal = 0;
    if (strategy === "mean")   fillVal = vals.reduce((a,b)=>a+b,0) / vals.length;
    if (strategy === "median") {
      const s = [...vals].sort((a,b)=>a-b);
      fillVal = s.length % 2 === 0 ? (s[s.length/2-1]+s[s.length/2])/2 : s[Math.floor(s.length/2)];
    }
    DATA.forEach(row => {
      const v = parseFloat(row[col]);
      if (isNaN(v) || row[col] === "" || row[col] === null || row[col] === undefined)
        row[col] = fillVal;
    });
  });
  renderDataTable();
  render();
}

function dropNullRows() {
  const numCols = COLS.filter(c => c.type === "num").map(c => c.name);
  const before  = DATA.length;
  DATA = DATA.filter(row =>
    numCols.every(col => { const v = parseFloat(row[col]); return !isNaN(v); })
  );
  renderDataTable();
  render();
  const dropped = before - DATA.length;
  if (dropped) alert(`Dropped ${dropped} row${dropped > 1 ? "s" : ""} with missing values.`);
}

// ── Render the data table ─────────────────────────────────────────────
function renderDataTable() {
  const container = document.getElementById("dataTableWrap");
  if (!container) return;

  const displayCols = COLS.slice(0, 20); // cap at 20 cols for readability
  const totalRows   = DATA.length;
  const totalPages  = Math.max(1, Math.ceil(totalRows / EXPLORER_PAGE_SIZE));
  _explorerPage     = Math.min(_explorerPage, totalPages - 1);
  const startIdx    = _explorerPage * EXPLORER_PAGE_SIZE;
  const pageRows    = DATA.slice(startIdx, startIdx + EXPLORER_PAGE_SIZE);

  // Build column summaries
  const summaries = displayCols.map(col => {
    if (col.type === "num") {
      const vals = DATA.map(r => parseFloat(r[col.name])).filter(v => !isNaN(v));
      const nulls = DATA.length - vals.length;
      const mean  = vals.length ? (vals.reduce((a,b)=>a+b,0)/vals.length).toFixed(2) : "—";
      return `min ${col.min?.toFixed(2) ?? "—"} · mean ${mean} · max ${col.max?.toFixed(2) ?? "—"}`
           + (nulls ? ` · <span style="color:#e74c3c">${nulls} nulls</span>` : "");
    } else {
      const vals  = DATA.map(r => r[col.name]).filter(v => v !== null && v !== "" && v !== undefined);
      const nulls = DATA.length - vals.length;
      const top   = [...new Set(vals)].slice(0,2).map(v => esc(String(v))).join(", ");
      return `${col.uniq} unique · ${top}${col.uniq > 2 ? "…" : ""}`
           + (nulls ? ` · <span style="color:#e74c3c">${nulls} nulls</span>` : "");
    }
  });

  const typeCycleMap = { cat:"num", num:"date", date:"cat" };
  const typeLabelMap = { cat:"cat", num:"num", date:"date" };

  const headerCells = displayCols.map((col, ci) => {
    const tb = col.type === "num" ? "tb-num" : (col.type === "date" ? "tb-date" : "tb-cat");
    return `<th>
      <div class="dex-col-header">
        <span class="dex-col-name" title="Double-click to rename"
              ondblclick="promptRenameCol('${esc(col.name)}')">${esc(col.name)}</span>
        <span class="type-badge ${tb}" style="cursor:pointer" title="Click to change type"
              onclick="overrideColType('${esc(col.name)}','${typeCycleMap[col.type]}')">${typeLabelMap[col.type]}</span>
      </div>
      <div class="dex-col-summary">${summaries[ci]}</div>
    </th>`;
  }).join("");

  const bodyRows = pageRows.map((row, ri) => {
    const cells = displayCols.map(col => {
      const v = row[col.name];
      const isEmpty = v === null || v === "" || v === undefined || (col.type === "num" && isNaN(parseFloat(v)));
      return `<td class="${isEmpty ? "dex-null" : ""}">${esc(String(v ?? ""))}</td>`;
    }).join("");
    return `<tr><td class="dex-rownum">${startIdx + ri + 1}</td>${cells}</tr>`;
  }).join("");

  // Pagination controls
  const pageInfo = `Showing rows ${startIdx+1}–${Math.min(startIdx+EXPLORER_PAGE_SIZE, totalRows)} of ${totalRows.toLocaleString()}`;
  const prevDis  = _explorerPage === 0 ? "disabled" : "";
  const nextDis  = _explorerPage >= totalPages-1 ? "disabled" : "";

  container.innerHTML = `
    <div class="dex-toolbar">
      <input type="text" class="dex-search" placeholder="🔍 Search rows…" value="${esc(ROW_FILTER_TEXT)}"
             oninput="applyRowFilter(this.value)">
      <div class="dex-actions">
        <span class="dex-action-lbl">Fill nulls:</span>
        <button class="dex-action-btn" onclick="fillNulls('zero')">→ 0</button>
        <button class="dex-action-btn" onclick="fillNulls('mean')">→ mean</button>
        <button class="dex-action-btn" onclick="fillNulls('median')">→ median</button>
        <button class="dex-action-btn dex-danger" onclick="dropNullRows()">Drop null rows</button>
      </div>
    </div>
    <div class="dex-table-wrap">
      <table class="dex-table">
        <thead><tr><th class="dex-rownum">#</th>${headerCells}</tr></thead>
        <tbody>${bodyRows}</tbody>
      </table>
    </div>
    <div class="dex-footer">
      <span class="dex-page-info">${pageInfo}${COLS.length > 20 ? ` (first 20 of ${COLS.length} columns shown)` : ""}</span>
      <div class="dex-pagination">
        <button class="dex-page-btn" onclick="explorerPage(-1)" ${prevDis}>← Prev</button>
        <span class="dex-page-num">Page ${_explorerPage+1} / ${totalPages}</span>
        <button class="dex-page-btn" onclick="explorerPage(1)" ${nextDis}>Next →</button>
      </div>
    </div>`;
}

function explorerPage(dir) {
  _explorerPage = Math.max(0, _explorerPage + dir);
  renderDataTable();
}

function promptRenameCol(oldName) {
  const newName = prompt(`Rename column "${oldName}" to:`, oldName);
  if (newName) renameCol(oldName, newName);
}
