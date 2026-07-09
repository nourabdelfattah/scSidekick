/*! © 2024–2026 Nourhan Abdelfattah — scSidekick R package
 *  Bundled for use via scSidekick::ChartBuilder() only.
 *  Not licensed for standalone redistribution.
 */
// ═══════════════════════════════════════════════════════════════════
//  CONSTANTS
// ═══════════════════════════════════════════════════════════════════
const FONT     = "-apple-system,'Helvetica Neue',Arial,sans-serif";
const MAX_ROWS = 10000;

const CHART_DEFS = [
  { id:"bar",       icon:"📊", name:"Bar",       req:{x:["cat","date"], y:["num"]} },
  { id:"scatter",   icon:"●",  name:"Scatter",   req:{x:["num"],        y:["num"]} },
  { id:"line",      icon:"📈", name:"Line",      req:{x:["any"],        y:["num"]} },
  { id:"histogram", icon:"▬",  name:"Histogram", req:{x:["num"],        y:[]} },
  { id:"box",       icon:"▭",  name:"Box",       req:{x:["cat"],        y:["num"]} },
  { id:"violin",    icon:"🎻", name:"Violin",    req:{x:["cat"],        y:["num"]} },
  { id:"heatmap",   icon:"▦",  name:"Heatmap",   req:{x:["cat"],        y:["cat"]} },
  { id:"alluvial",  icon:"≋",  name:"Alluvial",  req:{x:["cat"],        y:["cat"]} },
];

const CAT_PALS = {
  tableau: { name:"Tableau",     cb:false,
    colors:["#4e79a7","#f28e2b","#59a14f","#e15759","#76b7b2","#edc948","#b07aa1","#ff9da7"],
    stops: ["#4e79a7","#f28e2b","#59a14f","#e15759"] },
  okabe:   { name:"Okabe-Ito",   cb:true,
    colors:["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"],
    stops: ["#E69F00","#56B4E9","#009E73","#D55E00"] },
  nature:  { name:"Nature",      cb:false,
    colors:["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000"],
    stops: ["#E64B35","#4DBBD5","#00A087","#3C5488"] },
  d3:      { name:"D3 Cat10",    cb:false,
    colors:["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f"],
    stops: ["#1f77b4","#ff7f0e","#2ca02c","#d62728"] },
  bold:    { name:"Bold",        cb:false,
    colors:["#7b2d8b","#e63946","#f4a261","#2a9d8f","#264653","#e9c46a","#457b9d","#a8dadc"],
    stops: ["#7b2d8b","#e63946","#f4a261","#2a9d8f"] },
  pastel:  { name:"Pastel",      cb:false,
    colors:["#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec"],
    stops: ["#fbb4ae","#b3cde3","#ccebc5","#decbe4"] },
  mono:    { name:"Monochrome",  cb:false,
    colors:["#1a1a2e","#16213e","#0f3460","#533483","#2d6a4f","#1b4332","#6b4226","#7d5a50"],
    stops: ["#1a1a2e","#0f3460","#533483","#2d6a4f"] },
  viridis_cat: { name:"Viridis (cat)", cb:true,
    colors:["#440154","#443983","#31688e","#21908c","#35b779","#5dc963","#aadc32","#fde725"],
    stops: ["#440154","#31688e","#35b779","#fde725"] },
  plasma_cat:  { name:"Plasma (cat)",  cb:true,
    colors:["#0d0887","#5c01a6","#9c179e","#cc4778","#e16462","#f1924b","#fccd25","#f0f921"],
    stops: ["#0d0887","#9c179e","#e16462","#f0f921"] },
  custom:  { name:"Custom…",     cb:false,
    colors:[], stops:[] },
};

const SEQ_PALS = {
  viridis: { name:"Viridis",  cb:true,  scale:"Viridis",
    stops:["#440154","#31688e","#35b779","#fde725"] },
  plasma:  { name:"Plasma",   cb:true,  scale:"Plasma",
    stops:["#0d0887","#7e03a8","#f89441","#f0f921"] },
  magma:   { name:"Magma",    cb:true,  scale:"Magma",
    stops:["#000004","#51127c","#b73779","#fcfdbf"] },
  inferno: { name:"Inferno",  cb:true,  scale:"Inferno",
    stops:["#000004","#420a68","#d24b34","#fcffa4"] },
  cividis: { name:"Cividis",  cb:true,  scale:"Cividis",
    stops:["#00204d","#535d6c","#adbc44","#fee838"] },
  blues:   { name:"Blues",    cb:false, scale:"Blues",
    stops:["#eff3ff","#9ecae1","#3182bd","#08306b"] },
  greens:  { name:"Greens",   cb:false, scale:"Greens",
    stops:["#f7fcf5","#74c476","#238b45","#00441b"] },
  reds:    { name:"Reds",     cb:false, scale:"Reds",
    stops:["#fff5f0","#fc9272","#ef3b2c","#67000d"] },
  oranges: { name:"Oranges",  cb:false, scale:"Oranges",
    stops:["#fff5eb","#fd8d3c","#e6550d","#7f2704"] },
  teal:    { name:"Teal",     cb:false,
    scale:[[0,"#edf8f4"],[.4,"#74c8a3"],[.7,"#2a9d8f"],[1,"#0d3b36"]],
    stops:["#edf8f4","#74c8a3","#2a9d8f","#0d3b36"] },
};

const DIV_PALS = {
  rbu:      { name:"Blue–Red",        cb:true,
    scale:[[0,"#2166ac"],[0.25,"#92c5de"],[0.5,"#f7f7f7"],[0.75,"#f4a582"],[1,"#b2182b"]],
    stops:["#2166ac","#92c5de","#f7f7f7","#f4a582","#b2182b"] },
  rdylbu:   { name:"Red–Yellow–Blue", cb:true,  scale:"RdYlBu",
    stops:["#d73027","#fdae61","#ffffbf","#74add1","#4575b4"] },
  spectral: { name:"Spectral",        cb:false, scale:"Spectral",
    stops:["#d53e4f","#f46d43","#ffffbf","#66c2a5","#3288bd"] },
  coolwarm: { name:"Cool–Warm",       cb:true,
    scale:[[0,"#3b4cc0"],[0.25,"#88a7e6"],[0.5,"#f0f0f0"],[0.75,"#f5a47f"],[1,"#b40426"]],
    stops:["#3b4cc0","#88a7e6","#f0f0f0","#f5a47f","#b40426"] },
  bor:      { name:"Blue–Orange",     cb:false,
    scale:[[0,"#005f96"],[0.25,"#74b9e1"],[0.5,"#f7f7f7"],[0.75,"#f4a261"],[1,"#e65100"]],
    stops:["#005f96","#74b9e1","#f7f7f7","#f4a261","#e65100"] },
  pug:      { name:"Purple–Green",    cb:true,
    scale:[[0,"#762a83"],[0.25,"#c2a5cf"],[0.5,"#f7f7f7"],[0.75,"#7fbf7b"],[1,"#1b7837"]],
    stops:["#762a83","#c2a5cf","#f7f7f7","#7fbf7b","#1b7837"] },
};

// ═══════════════════════════════════════════════════════════════════
//  STATE
// ═══════════════════════════════════════════════════════════════════
let DATA = [], COLS = [];
let RAW_DATA  = [];          // original rows, never mutated
let FORMAT = "long";         // "long" | "wide"
let _pendingFilename = "";

// Data explorer state
let COL_OVERRIDES    = {};   // { colName: "num"|"cat"|"date" }
let COL_RENAMES      = {};   // { oldName: newName }
let ROW_FILTER_TEXT  = "";   // search string applied to RAW_DATA → DATA
let _explorerPage    = 0;    // current page in data table
const EXPLORER_PAGE_SIZE = 100;

const WS = {           // Wide-format state
  metaCols:     [],
  sampleCols:   [],
  sampleGroups: {},    // { sampleName: groupName }
  groupOrder:   [],    // ordered unique group names
  delimiter:    ".",
  sampleMeta:   {},    // { sampleName: { field1: val1, ... } }
  metaFields:   [],    // ordered list of metadata field names
  rowFilter:    null,  // null=all, or Set<string> of allowed rowLabel values
  sampleFilter: null,  // null=all, or Set<string> of allowed sample col names
};

const ST = {
  chartType:   "bar",
  barMode:     "group",
  catPal:      "tableau",
  seqPal:      "viridis",
  divPal:      "rbu",
  // Which scale family a heatmap currently uses: "seq" | "div"
  // Decoupled from the Center transform — Center just auto-switches this to "div".
  heatScaleType: "seq",
  // Custom categorical palette — 8 user-picked colors (default = Tableau)
  customColors: ["#4e79a7","#f28e2b","#59a14f","#e15759","#76b7b2","#edc948","#b07aa1","#ff9da7"],
  revX: false, revY: false,
  grid: false, labels: false,   // pub theme: no grid by default
  suggested:   [],
  log2:        false,
  center:      false,
  colSplit:    true,
  // Colorbar
  cbPosition:  "right",   // "right" | "bottom"
  // Figure style
  fontSize:    12,
  lineThick:   3,
  exportW:     1200,
  exportH:     800,
  // Faceting
  nFacetCols:  3,
  freeY:       false,
  // Individual data points on bar/box/violin
  showPoints:  "none",    // "none" | "outliers" | "all"
  blackPoints: true,      // true = jitter/dot points always black
  // Error bars on bar charts (only meaningful when aggFn = "mean")
  errorBars:   "sem",     // "none" | "sem" | "sd"
  // Axis line style — pub theme defaults (black axis lines, no grid)
  axisColor:   "#000000",
  axisThick:   1.5,       // 0 = hidden, >0 = line width in px
  // Y-axis range
  yFromZero:   false,     // true = rangemode:"tozero"
  // Per-group color overrides { groupName: "#hexcolor" }
  groupColors: {},
  // Long-format column value filters
  colFilters:  {},         // { colName: Set<string> }
};
