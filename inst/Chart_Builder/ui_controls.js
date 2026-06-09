/*! © 2024 Nourhan Abdelfattah — scSidekick R package
 *  Bundled for use via scSidekick::ChartBuilder() only.
 *  Not licensed for standalone redistribution.
 */
// ═══════════════════════════════════════════════════════════════════
//  UI CONTROLS  (depends on: config.js, statistics.js)
// ═══════════════════════════════════════════════════════════════════

const TYPE_DOT = { num:"#3b82f6", cat:"#22c55e", date:"#f59e0b" };

// ── Tab switching (Problem 4) ────────────────────────────────────────
function switchTab(tab) {
  document.getElementById("chartPanel").style.display = tab === "chart" ? "flex" : "none";
  document.getElementById("dataPanel").style.display  = tab === "data"  ? "flex" : "none";
  document.getElementById("statsPanel").style.display = tab === "stats" ? "flex" : "none";
  document.querySelectorAll(".main-tab-btn").forEach(b =>
    b.classList.toggle("on", b.dataset.tab === tab)
  );
  if (tab === "data")  renderDataTable();
}

// ── Column pills ─────────────────────────────────────────────────────
function buildColPills() {
  const c = document.getElementById("colPills");
  c.innerHTML = "";
  const effectiveCols = (FORMAT === "wide" && ST.chartType !== "heatmap")
    ? getMeltedCols()
    : (FORMAT === "wide"
        ? [...WS.metaCols.map(n => ({name:n, type:"cat"})),
           ...WS.sampleCols.slice(0,6).map(n => ({name:n, type:"num"})),
           ...(WS.sampleCols.length > 6 ? [{name:`+${WS.sampleCols.length-6} more`, type:"num"}] : [])]
        : COLS.slice(0, 16));
  effectiveCols.forEach(col => {
    const p = document.createElement("div");
    p.className = "col-pill";
    p.innerHTML = `<div class="col-dot" style="background:${TYPE_DOT[col.type]||'#aaa'}"></div><span>${esc(col.name)}</span>`;
    c.appendChild(p);
  });
}

// ── Chart type grid ──────────────────────────────────────────────────
function buildChartGrid() {
  const g = document.getElementById("chartGrid");
  g.innerHTML = "";
  CHART_DEFS.forEach(ct => {
    const avail = isChartAvailable(ct.id);
    const b = document.createElement("div");
    b.className = "chart-btn"
      + (ct.id === ST.chartType ? " on" : "")
      + (ST.suggested.includes(ct.id) ? " suggested" : "")
      + (avail ? "" : " unavailable");
    b.id = "cb_" + ct.id;
    b.title = avail ? "" : getUnavailableReason(ct.id);
    b.innerHTML = `<span class="chart-icon">${ct.icon}</span><span class="chart-name">${ct.name}</span>`;
    b.onclick = avail ? () => setChartType(ct.id) : null;
    g.appendChild(b);
  });
}

// Determine if a chart type is available given current data and format
function isChartAvailable(chartId) {
  if (chartId === "alluvial") {
    if (FORMAT === "wide") return false;   // alluvial needs tidy rows
    return COLS.filter(c => c.type === "cat").length >= 2;
  }
  if (chartId === "violin") {
    if (FORMAT === "wide") return true;    // melted data has cat + num
    return COLS.some(c => c.type === "cat") && COLS.some(c => c.type === "num");
  }
  if (FORMAT !== "wide") return true;      // long format: all other types available
  if (chartId === "heatmap") return true;
  if (chartId === "scatter") {
    const numMetaCols = COLS.filter(c => WS.metaCols.includes(c.name) && c.type === "num");
    return numMetaCols.length >= 1;
  }
  return true;
}
function getUnavailableReason(chartId) {
  if (chartId === "alluvial") {
    if (FORMAT === "wide") return "Alluvial needs long / tidy data. Switch format to Long first.";
    return "Alluvial needs at least 2 categorical columns (Source and Target).";
  }
  if (chartId === "scatter" && FORMAT === "wide")
    return "Scatter requires two numeric dimensions. No numeric annotation column detected.";
  return "";
}

// ── Palette selector (dropdown + gradient preview + optional custom picker) ──
// Two contexts only:
//   • heat — any heatmap (wide OR long). Dropdown offers ALL sequential + diverging
//            palettes, grouped. Center just auto-switches ST.heatScaleType to "div".
//   • cat  — every other chart type. Dropdown offers the categorical palettes (+ Custom).
function _palContext() {
  return { isHeat: ST.chartType === "heatmap" };
}

// Stops for any heat palette key (looks in both seq + div tables)
function _heatStops(key) { return (SEQ_PALS[key] || DIV_PALS[key])?.stops || []; }
function _heatActiveKey() { return ST.heatScaleType === "div" ? ST.divPal : ST.seqPal; }

function buildPalGrid() {
  const sel = document.getElementById("palSelect");
  if (!sel) return;
  const { isHeat } = _palContext();

  if (isHeat) {
    const active = _heatActiveKey();
    const optHtml = (dict) => Object.entries(dict).map(([key, p]) =>
      `<option value="${key}"${key === active ? " selected" : ""}>${p.name}${p.cb ? "  ✓CB" : ""}</option>`
    ).join("");
    sel.innerHTML =
      `<optgroup label="Sequential">${optHtml(SEQ_PALS)}</optgroup>` +
      `<optgroup label="Diverging">${optHtml(DIV_PALS)}</optgroup>`;
    _updatePalPreview(active);
    const customPicker = document.getElementById("customPalPicker");
    if (customPicker) customPicker.style.display = "none";
    return;
  }

  // Categorical context
  const active = ST.catPal;
  sel.innerHTML = Object.entries(CAT_PALS).map(([key, p]) =>
    `<option value="${key}"${key === active ? " selected" : ""}>${p.name}${p.cb ? "  ✓CB" : ""}</option>`
  ).join("");
  _updatePalPreview(active);

  const customPicker = document.getElementById("customPalPicker");
  if (customPicker) {
    const show = active === "custom";
    customPicker.style.display = show ? "" : "none";
    if (show) _buildCustomColorPicker();
  }
}

function onPalChange() {
  const { isHeat } = _palContext();
  const key = document.getElementById("palSelect").value;

  if (isHeat) {
    // Detect which family the chosen palette belongs to (keys never collide)
    if (DIV_PALS[key]) { ST.divPal = key; ST.heatScaleType = "div"; }
    else               { ST.seqPal = key; ST.heatScaleType = "seq"; }
    _updatePalPreview(key);
    render();
    return;
  }

  ST.catPal = key;
  _updatePalPreview(key);
  const customPicker = document.getElementById("customPalPicker");
  if (customPicker) {
    const show = key === "custom";
    customPicker.style.display = show ? "" : "none";
    if (show) _buildCustomColorPicker();
  }
  render();
}

function _updatePalPreview(key) {
  const strip = document.getElementById("palPreviewStrip");
  if (!strip) return;
  const { isHeat } = _palContext();
  let stops;
  if (isHeat)                     stops = _heatStops(key);
  else if (key === "custom")      stops = ST.customColors.slice(0, 4);
  else                            stops = CAT_PALS[key]?.stops || [];
  strip.style.background = stops.length
    ? `linear-gradient(to right,${stops.join(",")})`
    : "#e3e6ea";
}

function _buildCustomColorPicker() {
  const grid = document.getElementById("customColorGrid");
  if (!grid) return;
  grid.innerHTML = "";
  ST.customColors.forEach((color, i) => {
    const inp = document.createElement("input");
    inp.type  = "color";
    inp.value = color;
    inp.title = `Color ${i + 1}`;
    inp.oninput = e => {
      ST.customColors[i] = e.target.value;
      // Keep preview strip in sync
      const strip = document.getElementById("palPreviewStrip");
      if (strip) strip.style.background = `linear-gradient(to right,${ST.customColors.slice(0,4).join(",")})`;
      render();
    };
    grid.appendChild(inp);
  });
}

// ── Column selectors ──────────────────────────────────────────────────
function populateSelectors() {
  // For wide + non-heatmap, use melted cols; otherwise use COLS
  const effectiveCols = (FORMAT === "wide" && ST.chartType !== "heatmap")
    ? getMeltedCols()
    : COLS;

  const allOpts = effectiveCols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join("");
  const noneOpt = `<option value="">— None —</option>`;
  const catOpts = effectiveCols.filter(c => c.type === "cat")
    .map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join("");

  document.getElementById("selX").innerHTML     = allOpts;
  document.getElementById("selY").innerHTML     = allOpts;
  document.getElementById("selColor").innerHTML = noneOpt + allOpts;
  document.getElementById("selSize").innerHTML  = noneOpt + allOpts;
  document.getElementById("selFacet").innerHTML = noneOpt + catOpts;

  // Wide heatmap row selectors
  if (FORMAT === "wide") {
    const metaOpts = WS.metaCols
      .map(c => `<option value="${esc(c)}">${esc(c)}</option>`).join("");
    const rowLabelSel = document.getElementById("selRowLabel");
    const rowSplitSel = document.getElementById("selRowSplit");
    if (rowLabelSel) rowLabelSel.innerHTML = metaOpts || `<option value="">Row index</option>`;
    if (rowSplitSel) rowSplitSel.innerHTML = noneOpt + metaOpts;
    if (WS.metaCols.length >= 1 && rowLabelSel) rowLabelSel.value = WS.metaCols[0];
    // Problem 3.4: only auto-set row split when there are ≥2 meta cols
    if (WS.metaCols.length >= 2 && rowSplitSel)
      rowSplitSel.value = WS.metaCols[WS.metaCols.length - 1];
  }
}

function esc(s) {
  return String(s).replace(/&/g,"&amp;").replace(/"/g,"&quot;").replace(/</g,"&lt;");
}
// Safe embedding of a string inside a JS single-quoted string in an HTML onclick attribute
function jss(s) { return String(s).replace(/\\/g,"\\\\").replace(/'/g,"\\'"); }

// ═══════════════════════════════════════════════════════════════════
//  SMART DEFAULTS
// ═══════════════════════════════════════════════════════════════════
function applyDefaults() {
  const ct = ST.chartType;

  // Problem 5: smart defaults for wide + non-heatmap (melted data)
  if (FORMAT === "wide" && ct !== "heatmap") {
    const mc   = getMeltedCols();
    const grp  = mc.find(c => c.name === "Group");
    const smp  = mc.find(c => c.name === "Sample");
    const val  = mc.find(c => c.name === "Value");
    const meta = mc.find(c => WS.metaCols.includes(c.name));
    const xDefault = grp?.name || smp?.name;

    if (ct === "bar" || ct === "box" || ct === "line") {
      setVal(document.getElementById("selX"),     xDefault);
      setVal(document.getElementById("selY"),     val?.name);
      setVal(document.getElementById("selColor"), meta?.name || "");
    } else if (ct === "histogram") {
      setVal(document.getElementById("selX"),     val?.name);
      setVal(document.getElementById("selColor"), grp?.name || "");
    } else if (ct === "scatter") {
      const numMeta = mc.find(c => WS.metaCols.includes(c.name) && c.type === "num");
      setVal(document.getElementById("selX"),     numMeta?.name || val?.name);
      setVal(document.getElementById("selY"),     val?.name);
      setVal(document.getElementById("selColor"), grp?.name || "");
    }
    updateAxisUI();
    return;
  }

  // Long format defaults
  const cats  = COLS.filter(c => c.type === "cat");
  const nums  = COLS.filter(c => c.type === "num");
  const dates = COLS.filter(c => c.type === "date");

  const bestCat  = [...cats].sort((a,b) => a.uniq - b.uniq).find(c => c.uniq > 1) || cats[0];
  const bestCat2 = cats.find(c => c !== bestCat);
  const bestNum  = nums[0];
  const bestNum2 = nums[1];
  const bestDate = dates[0];

  const sx = document.getElementById("selX");
  const sy = document.getElementById("selY");
  const sc = document.getElementById("selColor");

  if (ct === "bar" || ct === "box" || ct === "violin") {
    setVal(sx, bestCat?.name); setVal(sy, bestNum?.name); setVal(sc, bestCat2?.name || "");
  } else if (ct === "scatter") {
    setVal(sx, bestNum?.name); setVal(sy, bestNum2?.name || bestNum?.name); setVal(sc, bestCat?.name || "");
  } else if (ct === "histogram") {
    setVal(sx, bestNum?.name); setVal(sc, bestCat?.name || "");
  } else if (ct === "heatmap") {
    setVal(sx, bestCat?.name); setVal(sy, bestCat2?.name || cats[0]?.name); setVal(sc, bestNum?.name || "");
  } else if (ct === "line") {
    setVal(sx, bestDate?.name || bestCat?.name || bestNum?.name);
    setVal(sy, bestNum?.name); setVal(sc, bestCat?.name || "");
  } else if (ct === "alluvial") {
    // Source = lowest-cardinality cat, Target = second cat, no middle node by default
    setVal(sx, bestCat?.name); setVal(sy, bestCat2?.name || ""); setVal(sc, "");
  }
  updateAxisUI();
}

function setVal(el, val) {
  if (!el || val === undefined || val === null) return;
  const opt = [...el.options].find(o => o.value === String(val));
  if (opt) el.value = String(val);
}

// ═══════════════════════════════════════════════════════════════════
//  UI STATE MANAGEMENT
// ═══════════════════════════════════════════════════════════════════
function updateAxisUI() {
  const ct = ST.chartType;
  const isWideHeat = FORMAT === "wide" && ct === "heatmap";

  document.getElementById("colMappingDivider").style.display = isWideHeat ? "none" : "";
  document.getElementById("colMappingSection").style.display = isWideHeat ? "none" : "";
  document.getElementById("wideHeatDivider").style.display   = isWideHeat ? "" : "none";
  document.getElementById("wideHeatSection").style.display   = isWideHeat ? "" : "none";
  document.getElementById("axisOptDivider").style.display    = isWideHeat ? "none" : "";
  document.getElementById("axisOptSection").style.display    = isWideHeat ? "none" : "";
  hide("colGroupOrderRow", !ST.colSplit || WS.groupOrder.length <= 1);

  // Filter panels
  const isWideFormat = FORMAT === "wide";
  hide("wideFilterDivider", !isWideFormat);
  hide("wideFilterSection",  !isWideFormat);
  hide("longFilterDivider",  isWideFormat);
  hide("longFilterSection",  isWideFormat);
  // Rebuild filter lists when needed
  if (isWideFormat) {
    if (typeof buildWideFilterPanel === "function") buildWideFilterPanel();
  } else {
    if (typeof buildLongFilterUI === "function") buildLongFilterUI();
  }

  if (isWideHeat) return;

  const hasColor  = !!document.getElementById("selColor").value;
  const hasFacet  = !!document.getElementById("selFacet")?.value;

  hide("yRow",           ct === "histogram");
  hide("sizeRow",        ct !== "scatter" && ct !== "alluvial");
  hide("binsRow",        ct !== "histogram");
  hide("showPointsRow",  !["bar","box","violin"].includes(ct));
  hide("blackPointsRow", !["bar","box","violin"].includes(ct) || ST.showPoints === "none");
  hide("errorBarRow",    ct !== "bar");
  hide("aggRow",         !["bar","line","box","heatmap"].includes(ct));
  hide("barModeRow",     ct !== "bar" || !hasColor);
  hide("sortRow",        ["scatter","histogram","heatmap","alluvial","violin"].includes(ct));
  // Alluvial can't be meaningfully faceted — hide those controls
  hide("facetRow",       ct === "alluvial");
  hide("facetNColsRow",  !hasFacet || ct === "alluvial");
  hide("freeYRow",       !hasFacet || ct === "alluvial");
  document.getElementById("swapBtn").style.display = ct === "heatmap" ? "none" : "";

  // Size-row label changes meaning for alluvial
  const sizeLbl = document.getElementById("sizeLbl");
  if (sizeLbl) sizeLbl.textContent = ct === "alluvial" ? "Value / weight" : "Bubble size";

  if (!COLS.length && !(FORMAT === "wide")) return;

  const effectiveCols = (FORMAT === "wide" && ct !== "heatmap") ? getMeltedCols() : COLS;
  const selColor  = document.getElementById("selColor");
  const curColor  = selColor.value;

  if (ct === "heatmap") {
    document.getElementById("colorLbl").textContent = "Z value (color scale)";
    const numOpts = effectiveCols.filter(c => c.type === "num")
      .map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join("");
    selColor.innerHTML = `<option value="">— None —</option>` + numOpts;
  } else if (ct === "alluvial") {
    document.getElementById("colorLbl").textContent = "Flow through / Middle node (optional)";
    const catOpts2 = effectiveCols.filter(c => c.type === "cat")
      .map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join("");
    selColor.innerHTML = `<option value="">— None (2-layer) —</option>` + catOpts2;
  } else {
    document.getElementById("colorLbl").textContent = "Color / Group by";
    const allOpts = effectiveCols.map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join("");
    selColor.innerHTML = `<option value="">— None —</option>` + allOpts;
  }
  setVal(selColor, curColor);

  const selFacet  = document.getElementById("selFacet");
  const curFacet  = selFacet?.value;
  const catOpts   = effectiveCols.filter(c => c.type === "cat")
    .map(c => `<option value="${esc(c.name)}">${esc(c.name)}</option>`).join("");
  if (selFacet) { selFacet.innerHTML = `<option value="">— None —</option>` + catOpts; setVal(selFacet, curFacet); }

  const aggLbl = document.querySelector("#aggRow .ctrl-lbl");
  if (aggLbl) aggLbl.textContent = ct === "heatmap" ? "Aggregate Z by" : "Aggregate Y by";

  const xHint = {bar:"Category",scatter:"Number",line:"Any",histogram:"Number",box:"Category",violin:"Category",heatmap:"Category",alluvial:"Source"}[ct]||"Any";
  const yHint = {bar:"Number",scatter:"Number",line:"Number",histogram:"—",box:"Number",violin:"Number",heatmap:"Category",alluvial:"Target"}[ct]||"Number";
  const xCls  = {bar:"tb-cat",scatter:"tb-num",line:"tb-any",histogram:"tb-num",box:"tb-cat",violin:"tb-cat",heatmap:"tb-cat",alluvial:"tb-cat"}[ct]||"tb-any";
  const yCls  = {bar:"tb-num",scatter:"tb-num",line:"tb-num",histogram:"",box:"tb-num",violin:"tb-num",heatmap:"tb-cat",alluvial:"tb-cat"}[ct]||"tb-num";

  document.getElementById("xBadge").textContent = xHint;
  document.getElementById("xBadge").className   = "type-badge " + xCls;
  document.getElementById("yBadge").textContent = yHint;
  document.getElementById("yBadge").className   = "type-badge " + yCls;

  // Show "Run Stats" button for chart types that support statistics
  const statsEligible = ["bar","box","violin","scatter"].includes(ct);
  const statsBtn = document.getElementById("runStatsBtn");
  if (statsBtn) statsBtn.style.display = statsEligible ? "" : "none";

  // Per-group color swatches (only when a color column is selected)
  buildGroupColorSwatches();
}

function hide(id, condition) {
  const el = document.getElementById(id);
  if (el) el.style.display = condition ? "none" : "";
}

function setChartType(type) {
  ST.chartType = type;
  document.querySelectorAll(".chart-btn").forEach(b => b.classList.remove("on"));
  const btn = document.getElementById("cb_" + type);
  if (btn) btn.classList.add("on");
  buildPalGrid();
  // Re-populate selectors when switching between heatmap and other types in wide mode
  if (FORMAT === "wide") { populateSelectors(); buildColPills(); }
  applyDefaults();
  // Clear stats entirely if switching to a non-stats chart type;
  // otherwise show the rerun banner so user can decide
  const statsEligible = ["bar","box","violin","scatter"].includes(type);
  if (typeof SW !== "undefined" && SW.applied) {
    if (!statsEligible) {
      if (typeof clearStats === "function") clearStats();
    } else {
      if (typeof onAxisChangeWithStats === "function") onAxisChangeWithStats();
    }
  }
  render();
}

function setBarMode(mode, el) {
  ST.barMode = mode;
  document.querySelectorAll("#barModeRow .pill").forEach(p => p.classList.remove("on"));
  el.classList.add("on");
  render();
}

function onColorChange() {
  const hasColor = !!document.getElementById("selColor").value;
  hide("barModeRow", ST.chartType !== "bar" || !hasColor);
  if (typeof onAxisChangeWithStats === "function") onAxisChangeWithStats();
  render();
}

function onFacetChange() {
  const hasFacet = !!document.getElementById("selFacet")?.value;
  hide("facetNColsRow", !hasFacet);
  hide("freeYRow",      !hasFacet);
  render();
}

function onXYChange() {
  if (typeof onAxisChangeWithStats === "function") onAxisChangeWithStats();
  render();
}

function toggleSw(el, key) {
  el.classList.toggle("on");
  ST[key] = el.classList.contains("on");
  if (key === "center") {
    // Centering makes data symmetric around 0 → auto-switch to a diverging scale.
    // (User can still override the palette afterward via the dropdown.)
    ST.heatScaleType = ST.center ? "div" : "seq";
    buildPalGrid();
  }
  render();
}

function toggleColSplit() {
  const sw = document.getElementById("togColSplit");
  sw.classList.toggle("on");
  ST.colSplit = sw.classList.contains("on");
  hide("colGroupOrderRow", !ST.colSplit || WS.groupOrder.length <= 1);
  render();
}

function setShowPoints(mode) {
  ST.showPoints = mode;
  document.querySelectorAll("#showPointsRow .pill").forEach(p =>
    p.classList.toggle("on", p.dataset.pts === mode));
  render();
}

function setErrorBars(mode) {
  ST.errorBars = mode;
  document.querySelectorAll("#errorBarRow .pill").forEach(p =>
    p.classList.toggle("on", p.dataset.eb === mode));
  render();
}

function setBlackPoints(v) {
  ST.blackPoints = v;
  document.querySelectorAll("#blackPointsRow .pill").forEach(p =>
    p.classList.toggle("on", p.dataset.bp === (v ? "black" : "group")));
  render();
}

function setYFromZero(el) {
  el.classList.toggle("on");
  ST.yFromZero = el.classList.contains("on");
  render();
}

function setGroupColor(name, hex) {
  if (!ST.groupColors) ST.groupColors = {};
  ST.groupColors[name] = hex;
  render();
}

function buildGroupColorSwatches() {
  const section = document.getElementById("groupColorsSection");
  const list    = document.getElementById("groupColorsList");
  if (!section || !list) return;

  const colorCol = document.getElementById("selColor")?.value || null;
  if (!colorCol || !DATA.length) { section.style.display = "none"; return; }

  // Determine groups: wide format "Group" column → WS.groupOrder, else read from DATA
  let groups;
  if (FORMAT === "wide" && colorCol === "Group" && WS.groupOrder.length) {
    groups = WS.groupOrder;
  } else if (DATA[0] && colorCol in DATA[0]) {
    groups = [...new Set(DATA.map(d => String(d[colorCol] ?? "")))].filter(Boolean);
  } else {
    // colorCol is a melted pseudo-column; not available directly in DATA
    section.style.display = "none"; return;
  }

  if (!groups.length || groups.length > 30) { section.style.display = "none"; return; }

  const colors = getColors();
  section.style.display = "";
  list.innerHTML = groups.map((g, gi) => {
    const current = (ST.groupColors && ST.groupColors[g]) || colors[gi % colors.length];
    return `<div style="display:flex;align-items:center;gap:7px;padding:1px 0">
      <input type="color" value="${current}"
             oninput="setGroupColor('${jss(g)}', this.value)"
             style="width:22px;height:18px;border:none;border-radius:4px;cursor:pointer;padding:0;background:none;flex-shrink:0">
      <span style="font-size:11px;color:#3d4a5c;overflow:hidden;text-overflow:ellipsis;white-space:nowrap">${esc(g)}</span>
    </div>`;
  }).join("");
}

function setAxisColor(val) {
  ST.axisColor = val;
  render();
}

function setAxisThick(val) {
  ST.axisThick = parseFloat(val);
  document.getElementById("axisThickVal").textContent = val + "px";
  render();
}

function swapAxes() {
  const sx = document.getElementById("selX");
  const sy = document.getElementById("selY");
  const tmp = sx.value; sx.value = sy.value; sy.value = tmp;
  render();
}

// Problem 6: colorbar position toggle
// Must purge + fresh-draw because Plotly.react() doesn't reliably restructure
// the colorbar SVG when switching between vertical (right) and horizontal (bottom).
function setCbPosition(pos) {
  ST.cbPosition = pos;
  document.querySelectorAll(".cb-pos-btn").forEach(b =>
    b.classList.toggle("on", b.dataset.pos === pos)
  );
  Plotly.purge("plot");
  requestAnimationFrame(() => render());
}

// Problem 7: figure style controls
function setFontSize(val) {
  ST.fontSize = parseInt(val);
  document.getElementById("fontSizeVal").textContent = val + "px";
  render();
}
function setLineThick(val) {
  ST.lineThick = parseFloat(val);
  document.getElementById("lineThickVal").textContent = val + "px";
  render();
}
function setExportPreset(w, h) {
  document.getElementById("exportW").value = w;
  document.getElementById("exportH").value = h;
  ST.exportW = w; ST.exportH = h;
}
function onExportDimChange() {
  ST.exportW = parseInt(document.getElementById("exportW").value) || 2400;
  ST.exportH = parseInt(document.getElementById("exportH").value) || 1600;
}
function toggleStyleSection() {
  const body     = document.getElementById("styleSectionBody");
  const ico      = document.getElementById("styleToggleIco");
  // Check computed display (not inline style) so the CSS default of display:none is seen
  const isHidden = getComputedStyle(body).display === "none";
  body.style.display = isHidden ? "block" : "none";
  ico.textContent    = isHidden ? "▾" : "▸";
}

// ═══════════════════════════════════════════════════════════════════
//  MISC
// ═══════════════════════════════════════════════════════════════════
function showSpin(msg) {
  document.getElementById("spinMsg").textContent = msg || "Processing…";
  document.getElementById("spinner").classList.add("show");
}
function hideSpin() { document.getElementById("spinner").classList.remove("show"); }

// ═══════════════════════════════════════════════════════════════════
//  FEATURE 3 — DATA FILTER PANEL
// ═══════════════════════════════════════════════════════════════════

// ─── Wide format: Gene / row filter ──────────────────────────────────────────
function buildWideFilterPanel() {
  updateGeneFilterList(document.getElementById("geneSearchInput")?.value || "");
  buildSampleFilterList();
}

function updateGeneFilterList(search) {
  search = search || "";
  const rowLabelCol = document.getElementById("selRowLabel")?.value || WS.metaCols[0] || null;
  if (!rowLabelCol || !DATA.length) return;

  const allGenes = [...new Set(DATA.map(r => String(r[rowLabelCol] ?? "")))].filter(Boolean).sort();
  const filtered  = search
    ? allGenes.filter(g => g.toLowerCase().includes(search.toLowerCase()))
    : allGenes;

  const selected = WS.rowFilter || new Set(allGenes);
  const list     = document.getElementById("geneFilterList");
  if (!list) return;

  const shown = filtered.slice(0, 300);
  list.innerHTML = shown.map(g =>
    `<label class="filter-item">
      <input type="checkbox" ${selected.has(g) ? "checked" : ""}
             onchange="_toggleGene('${esc(g)}',this.checked)">
      <span>${esc(g)}</span>
    </label>`
  ).join("")
    + (filtered.length > 300
      ? `<div style="font-size:11px;color:#9aa5b4;padding:3px 6px">…${filtered.length - 300} more — refine search</div>`
      : "");

  const selCount = WS.rowFilter ? WS.rowFilter.size : allGenes.length;
  const status = document.getElementById("geneFilterStatus");
  if (status) status.textContent = (WS.rowFilter && WS.rowFilter.size < allGenes.length)
    ? `${selCount} of ${allGenes.length} selected`
    : `All ${allGenes.length} selected`;
  _updateFilterActiveTag();
}

function _toggleGene(name, checked) {
  const rowLabelCol = document.getElementById("selRowLabel")?.value || WS.metaCols[0] || null;
  if (!rowLabelCol) return;
  const allGenes = new Set(DATA.map(r => String(r[rowLabelCol] ?? "")).filter(Boolean));
  if (!WS.rowFilter) WS.rowFilter = new Set(allGenes);
  if (checked) WS.rowFilter.add(name); else WS.rowFilter.delete(name);
  if (WS.rowFilter.size >= allGenes.size) WS.rowFilter = null;
  const selCount = WS.rowFilter ? WS.rowFilter.size : allGenes.size;
  const status = document.getElementById("geneFilterStatus");
  if (status) status.textContent = WS.rowFilter ? `${selCount} of ${allGenes.size} selected` : `All ${allGenes.size} selected`;
  _updateFilterActiveTag();
  _filterChanged();
}

function selectAllGenes() {
  WS.rowFilter = null;
  updateGeneFilterList(document.getElementById("geneSearchInput")?.value || "");
  _filterChanged();
}

function clearAllGenes() {
  WS.rowFilter = new Set();
  updateGeneFilterList(document.getElementById("geneSearchInput")?.value || "");
  _filterChanged();
}

function togglePasteGeneArea() {
  const el = document.getElementById("pasteGeneArea");
  if (el) el.style.display = el.style.display === "none" ? "" : "none";
}

function applyPastedGeneList() {
  const ta = document.getElementById("geneListPaste");
  if (!ta) return;
  const names = ta.value.split(/[\n,;]+/).map(s => s.trim()).filter(Boolean);
  if (!names.length) return;
  const rowLabelCol = document.getElementById("selRowLabel")?.value || WS.metaCols[0] || null;
  if (!rowLabelCol) return;
  const allGenes = new Set(DATA.map(r => String(r[rowLabelCol] ?? "")).filter(Boolean));
  const matched  = names.filter(n => allGenes.has(n));
  WS.rowFilter   = new Set(matched);
  const notFound = names.filter(n => !allGenes.has(n));
  const status = document.getElementById("geneFilterStatus");
  if (status && notFound.length) {
    status.textContent = `${matched.length} matched · ${notFound.length} not found: ${notFound.slice(0,3).join(", ")}${notFound.length>3?"…":""}`;
  }
  togglePasteGeneArea();
  updateGeneFilterList("");
  _filterChanged();
}

// ─── Wide format: Sample filter ───────────────────────────────────────────────
let _sampleFilterMode = "group";

function setSampleFilterMode(mode) {
  _sampleFilterMode = mode;
  document.getElementById("sampleFilterModeGroup")?.classList.toggle("on",  mode === "group");
  document.getElementById("sampleFilterModeSample")?.classList.toggle("on", mode === "sample");
  buildSampleFilterList();
}

function buildSampleFilterList() {
  const list = document.getElementById("sampleFilterList");
  if (!list) return;

  if (_sampleFilterMode === "group") {
    const groups   = WS.groupOrder.length ? WS.groupOrder : [...new Set(Object.values(WS.sampleGroups))];
    const selected = WS.sampleFilter;
    list.innerHTML = groups.map(g => {
      const samplesInGroup = WS.sampleCols.filter(sc => WS.sampleGroups[sc] === g);
      const isChecked = !selected || samplesInGroup.some(sc => selected.has(sc));
      return `<label class="filter-item">
        <input type="checkbox" ${isChecked ? "checked" : ""}
               onchange="_toggleSampleGroup('${esc(g)}',this.checked)">
        <span>${esc(g)} <small style="color:#9aa5b4">(${samplesInGroup.length})</small></span>
      </label>`;
    }).join("");
  } else {
    const selected = WS.sampleFilter;
    list.innerHTML = WS.sampleCols.map(sc => {
      const isChecked = !selected || selected.has(sc);
      return `<label class="filter-item">
        <input type="checkbox" ${isChecked ? "checked" : ""}
               onchange="_toggleSample('${esc(sc)}',this.checked)">
        <span style="font-size:11px">${esc(sc)}</span>
      </label>`;
    }).join("");
  }
}

function _toggleSampleGroup(groupName, checked) {
  const samplesInGroup = WS.sampleCols.filter(sc => WS.sampleGroups[sc] === groupName);
  if (!WS.sampleFilter) WS.sampleFilter = new Set(WS.sampleCols);
  if (checked) samplesInGroup.forEach(sc => WS.sampleFilter.add(sc));
  else         samplesInGroup.forEach(sc => WS.sampleFilter.delete(sc));
  if (WS.sampleFilter.size >= WS.sampleCols.length) WS.sampleFilter = null;
  _updateFilterActiveTag();
  _filterChanged();
}

function _toggleSample(sampleName, checked) {
  if (!WS.sampleFilter) WS.sampleFilter = new Set(WS.sampleCols);
  if (checked) WS.sampleFilter.add(sampleName);
  else         WS.sampleFilter.delete(sampleName);
  if (WS.sampleFilter.size >= WS.sampleCols.length) WS.sampleFilter = null;
  buildSampleFilterList();
  _updateFilterActiveTag();
  _filterChanged();
}

function selectAllSamples() {
  WS.sampleFilter = null;
  buildSampleFilterList();
  _updateFilterActiveTag();
  _filterChanged();
}

function clearAllSamples() {
  WS.sampleFilter = new Set();
  buildSampleFilterList();
  _updateFilterActiveTag();
  _filterChanged();
}

// ─── Long format: column value filters ────────────────────────────────────────
// Columns with more than this many unique values get the search-driven UI
const LONG_HIGH_CARD = 20;

function buildLongFilterUI() {
  const container = document.getElementById("longFilterContent");
  if (!container) return;
  // No cardinality cap — show ALL categorical columns
  const catCols = COLS.filter(c => c.type === "cat" && c.uniq > 0);
  if (!catCols.length) {
    container.innerHTML = `<div style="font-size:12px;color:#9aa5b4">No categorical columns to filter.</div>`;
    return;
  }
  container.innerHTML = catCols.map((col, ci) => {
    const vals     = [...new Set(DATA.map(r => String(r[col.name] ?? "")))].filter(Boolean).sort();
    const active   = ST.colFilters[col.name];
    const allSel   = !active;                       // null = all; empty Set = none
    const highCard = col.uniq > LONG_HIGH_CARD;
    const cn       = jss(col.name);
    const countTxt = allSel ? "All" : (active.size + "/" + vals.length);
    const countClr = allSel ? "#9aa5b4" : "#4a8fe8";

    return `<div>
      <div class="filter-col-head"
           onclick="var v=document.getElementById('fcv_${ci}');v.style.display=v.style.display==='none'?'':'none'">
        <span>${esc(col.name)}</span>
        <span id="fcc_${ci}" style="font-size:10px;color:${countClr}">${countTxt}</span>
      </div>
      <div class="filter-col-vals" id="fcv_${ci}" style="display:none">
        <div style="display:flex;gap:5px;margin-bottom:4px;flex-wrap:wrap">
          <button class="pill on" style="font-size:10px;padding:2px 7px" onclick="_filterColAll('${cn}',${ci})">All</button>
          <button class="pill"    style="font-size:10px;padding:2px 7px" onclick="_filterColNone('${cn}',${ci})">None</button>
          ${highCard ? `<button class="pill" style="font-size:10px;padding:2px 7px" onclick="_togglePasteLongArea(${ci})">Paste list</button>` : ""}
        </div>
        ${highCard ? `<input type="text" id="fcs_${ci}" class="ax-input"
          style="width:100%;padding:5px 8px;border:1px solid #dde1e7;border-radius:7px;font-size:12px;margin-bottom:4px;box-sizing:border-box"
          placeholder="Search ${esc(col.name)}…"
          oninput="_searchLongCol('${cn}',${ci},this.value)">` : ""}
        <div id="fcvlist_${ci}" style="max-height:160px;overflow-y:auto;display:flex;flex-direction:column;gap:2px">
          ${_buildLongCheckboxes(col.name, vals, active, allSel, "")}
        </div>
        ${highCard ? `
        <div id="fcvpaste_${ci}" style="display:none;margin-top:4px">
          <textarea id="fcvtextarea_${ci}" rows="4"
            style="width:100%;padding:6px 9px;border:1px solid #dde1e7;border-radius:7px;font-size:12px;resize:vertical;box-sizing:border-box"
            placeholder="Paste one value per line…"></textarea>
          <button class="btn btn-blue" style="margin-top:4px;font-size:11px"
                  onclick="_applyPasteLongList('${cn}',${ci})">Apply</button>
        </div>` : ""}
      </div>
    </div>`;
  }).join("<hr style='border:none;border-top:1px solid #edf0f4;margin:4px 0'>");
}

// Build the inner checkbox HTML for one column (used by buildLongFilterUI & _searchLongCol)
function _buildLongCheckboxes(colName, vals, active, allSel, search) {
  const filtered = search
    ? vals.filter(v => v.toLowerCase().includes(search.toLowerCase()))
    : vals;
  const shown = filtered.slice(0, 300);
  const cn    = jss(colName);
  return shown.map(v => `
    <label class="filter-item">
      <input type="checkbox" ${(allSel || active?.has(v)) ? "checked" : ""}
             onchange="_toggleColFilter('${cn}','${jss(v)}',this.checked)">
      <span>${esc(v)}</span>
    </label>`).join("")
    + (filtered.length > 300
      ? `<div style="font-size:11px;color:#9aa5b4;padding:3px 6px">…${filtered.length - 300} more — refine search</div>`
      : "");
}

// Rebuild only the checkbox list for one column (preserves the rest of the panel)
function _searchLongCol(colName, ci, search) {
  const vals   = [...new Set(DATA.map(r => String(r[colName] ?? "")))].filter(Boolean).sort();
  const active = ST.colFilters[colName];
  const allSel = !active;
  const listEl = document.getElementById(`fcvlist_${ci}`);
  if (listEl) listEl.innerHTML = _buildLongCheckboxes(colName, vals, active, allSel, search);
}

// Refresh just the count header + checkbox list for one column (no full rebuild)
function _refreshLongColUI(colName, ci) {
  const vals   = [...new Set(DATA.map(r => String(r[colName] ?? "")))].filter(Boolean);
  const active = ST.colFilters[colName];
  const allSel = !active;
  const countEl = document.getElementById(`fcc_${ci}`);
  if (countEl) {
    countEl.textContent = allSel ? "All" : (active.size + "/" + vals.length);
    countEl.style.color = allSel ? "#9aa5b4" : "#4a8fe8";
  }
  const searchEl = document.getElementById(`fcs_${ci}`);
  _searchLongCol(colName, ci, searchEl?.value || "");
}

function _toggleColFilter(colName, val, checked) {
  const vals = [...new Set(DATA.map(r => String(r[colName] ?? "")))].filter(Boolean);
  if (!ST.colFilters[colName]) ST.colFilters[colName] = new Set(vals);
  if (checked) ST.colFilters[colName].add(val); else ST.colFilters[colName].delete(val);
  if (ST.colFilters[colName].size >= vals.length) delete ST.colFilters[colName];
  const ci = COLS.filter(c => c.type === "cat" && c.uniq > 0).findIndex(c => c.name === colName);
  if (ci >= 0) _refreshLongColUI(colName, ci);
  _updateFilterActiveTag();
  _filterChanged();
}

function _filterColAll(colName, ci) {
  delete ST.colFilters[colName];
  if (ci != null) _refreshLongColUI(colName, ci);
  _updateFilterActiveTag(); _filterChanged();
}
function _filterColNone(colName, ci) {
  ST.colFilters[colName] = new Set();
  if (ci != null) _refreshLongColUI(colName, ci);
  _updateFilterActiveTag(); _filterChanged();
}

// Paste-list helpers for high-cardinality columns
function _togglePasteLongArea(ci) {
  const el = document.getElementById(`fcvpaste_${ci}`);
  if (el) el.style.display = el.style.display === "none" ? "" : "none";
}
function _applyPasteLongList(colName, ci) {
  const ta = document.getElementById(`fcvtextarea_${ci}`);
  if (!ta) return;
  const names = ta.value.split(/[\n,;]+/).map(s => s.trim()).filter(Boolean);
  if (!names.length) return;
  const allVals = new Set(DATA.map(r => String(r[colName] ?? "")).filter(Boolean));
  const matched  = names.filter(n => allVals.has(n));
  ST.colFilters[colName] = new Set(matched);
  if (!ST.colFilters[colName].size || ST.colFilters[colName].size >= allVals.size)
    delete ST.colFilters[colName];
  _togglePasteLongArea(ci);
  buildLongFilterUI();  // full rebuild after paste (search is cleared anyway)
  _updateFilterActiveTag(); _filterChanged();
}

function clearAllColFilters() { ST.colFilters = {}; buildLongFilterUI(); _updateFilterActiveTag(); _filterChanged(); }

// ─── Shared filter utilities ──────────────────────────────────────────────────
function _filterChanged() {
  if (typeof SW !== "undefined" && SW.applied) {
    const banner = document.getElementById("statsRerunBanner");
    if (banner) banner.style.display = "flex";
  }
  render();
}

function _updateFilterActiveTag() {
  // Wide
  const rowLabelCol = document.getElementById("selRowLabel")?.value || WS.metaCols[0] || null;
  const allGenesCount = rowLabelCol && DATA.length
    ? new Set(DATA.map(r => String(r[rowLabelCol] ?? "")).filter(Boolean)).size
    : Infinity;
  const wideActive = (WS.rowFilter != null && WS.rowFilter.size < allGenesCount) || WS.sampleFilter != null;
  const wideTag = document.getElementById("filterActiveTag");
  if (wideTag) wideTag.style.display = wideActive ? "" : "none";
  // Long
  const longActive = Object.keys(ST.colFilters || {}).length > 0;
  const longTag = document.getElementById("longFilterActiveTag");
  if (longTag) longTag.style.display = longActive ? "" : "none";
}

function resetApp() {
  DATA = []; COLS = []; RAW_DATA = [];
  FORMAT = "long";
  COL_OVERRIDES = {}; COL_RENAMES = {}; ROW_FILTER_TEXT = ""; _explorerPage = 0;
  WS.metaCols = []; WS.sampleCols = []; WS.sampleGroups = {};
  WS.groupOrder = []; WS.delimiter = ".";
  WS.sampleMeta = {}; WS.metaFields = [];
  WS.rowFilter = null; WS.sampleFilter = null;
  ST.colFilters = {};
  _pendingFilename = "";
  Plotly.purge("plot");
  document.getElementById("appScreen").classList.remove("visible");
  document.getElementById("setupScreen").classList.remove("visible");
  document.getElementById("uploadScreen").style.display = "";
  fileInput.value = "";
}
