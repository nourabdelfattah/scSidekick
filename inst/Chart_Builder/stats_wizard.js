/*! © 2024–2026 Nourhan Abdelfattah — scSidekick R package
 *  Bundled for use via scSidekick::ChartBuilder() only.
 *  Not licensed for standalone redistribution.
 */
// ═══════════════════════════════════════════════════════════════════
//  STATS WIZARD  (depends on: config.js, statistics.js)
// ═══════════════════════════════════════════════════════════════════

// ── State ────────────────────────────────────────────────────────────
const SW = {
  open:         false,
  step:         0,      // 0-based wizard step index
  // detected context
  chartType:    null,   // "twoGroup"|"multiGroup"|"scatter"
  groups:       [],     // [{ label, values[] }]
  xVals:        [],     // scatter X
  yVals:        [],     // scatter Y
  xCol:         "",
  yCol:         "",
  colorCol:     "",
  // user choices
  testName:     null,   // "welch"|"paired"|"mannwhitney"|"anova"|"kruskal"|"pearson"|"spearman"
  posthoc:      "bonferroni", // "tukey"|"bonferroni"|"bh"|"none"
  compStyle:    "pairwise",   // "pairwise"|"reference"
  refGroup:     null,
  paired:       false,
  tails:        2,
  pDisplay:     "exact",      // "exact"|"stars"
  // normality results per group
  normality:    [],
  // final results
  results:      null,
  // flag: stats have been applied to current column mapping
  applied:      false,
};

// ── Open / close ─────────────────────────────────────────────────────
function openStatsWizard() {
  _detectContext();
  SW.step = 0;
  SW.open = true;
  _renderWizardStep();
  document.getElementById("statsWizardOverlay").style.display = "flex";
}

function closeStatsWizard() {
  SW.open = false;
  document.getElementById("statsWizardOverlay").style.display = "none";
}

function clearStats() {
  SW.applied  = false;
  SW.results  = null;
  SW.groups   = [];
  _clearBrackets();
  _clearStatsTab();
  // Clear only the stats section; keep the figure description
  const statsTa = document.getElementById("legendStatsTextarea");
  if (statsTa) statsTa.value = "";
  closeStatsWizard();
}

// ── Column-change warning ─────────────────────────────────────────────
function onAxisChangeWithStats() {
  if (!SW.applied) return;
  document.getElementById("statsRerunBanner").style.display = "flex";
}

function dismissRerunBanner(action) {
  document.getElementById("statsRerunBanner").style.display = "none";
  if (action === "rerun")  { clearStats(); openStatsWizard(); }
  if (action === "clear")  { clearStats(); }
  // "keep" — do nothing, stale brackets remain
}

// ── Context detection ────────────────────────────────────────────────
function _detectContext() {
  const ct = ST.chartType;
  SW.xCol    = document.getElementById("selX")?.value     || "";
  SW.yCol    = document.getElementById("selY")?.value     || "";
  SW.colorCol= document.getElementById("selColor")?.value || "";

  // For wide format, use the already-melted long data so groups/values are real
  const data = (FORMAT === "wide" && ct !== "heatmap") ? wideToLong() : DATA;

  if (ct === "scatter") {
    SW.chartType = "scatter";
    SW.xVals = data.map(d=>parseFloat(d[SW.xCol])).filter(v=>isFinite(v));
    SW.yVals = data.map(d=>parseFloat(d[SW.yCol])).filter(v=>isFinite(v));
    return;
  }

  // Groups come from colorCol if set, otherwise from xCol
  const groupKey = SW.colorCol || SW.xCol;
  const labels   = [...new Set(data.map(d=>String(d[groupKey])))];

  SW.groups = labels
    .map(lbl => ({
      label:  lbl,
      values: data.filter(d=>String(d[groupKey])===lbl)
                  .map(d=>parseFloat(d[SW.yCol])).filter(v=>isFinite(v)),
    }))
    .filter(g => g.values.length > 0);

  SW.chartType = SW.groups.length === 2 ? "twoGroup" : "multiGroup";
  SW.refGroup  = SW.groups[0]?.label || null;

  // Run Shapiro-Wilk per group (silently)
  SW.normality = SW.groups.map(g => ({
    label: g.label,
    n:     g.values.length,
    ...shapiroWilk(g.values),
  }));
}

// ── Default test recommendation ──────────────────────────────────────
function _recommendTest() {
  const allNormal = SW.normality.every(r => r.p > 0.05 || r.n < 3);
  if (SW.chartType === "scatter")    return allNormal ? "pearson" : "spearman";
  if (SW.chartType === "twoGroup")   return allNormal ? "welch"   : "mannwhitney";
  if (SW.chartType === "multiGroup") return allNormal ? "anova"   : "kruskal";
  return "welch";
}

// ── Wizard step renderer ──────────────────────────────────────────────
function _renderWizardStep() {
  const body = document.getElementById("wizardBody");
  const prev = document.getElementById("wizardPrevBtn");
  const next = document.getElementById("wizardNextBtn");

  prev.style.display = SW.step === 0 ? "none" : "";
  next.textContent   = SW.step >= _totalSteps()-1 ? "Run ▶" : "Next →";

  body.innerHTML = "";

  switch(SW.step) {
    case 0: _stepContext(body);    break;
    case 1: _stepNormality(body); break;
    case 2: _stepTestSelect(body);break;
    case 3: _stepOptions(body);   break;
    case 4: _stepConfirm(body);   break;
  }
}

function _totalSteps() {
  // Scatter skips options step (step 3)
  return SW.chartType === "scatter" ? 4 : 5;
}

function wizardNext() {
  const last = SW.step >= _totalSteps()-1;
  if (last) { _runStats(); return; }
  // Skip step 3 (options) for scatter
  if (SW.step === 2 && SW.chartType === "scatter") SW.step++;
  SW.step++;
  _renderWizardStep();
}

function wizardPrev() {
  if (SW.step === 0) return;
  if (SW.step === 4 && SW.chartType === "scatter") SW.step--;
  SW.step--;
  _renderWizardStep();
}

// Step 0 — Context
function _stepContext(body) {
  const ct = SW.chartType;
  const desc = ct === "scatter"
    ? `Correlation between <b>${SW.xCol}</b> and <b>${SW.yCol}</b>`
    : `Comparing <b>${SW.yCol}</b> across <b>${SW.groups.length} groups</b> on <b>${SW.xCol}</b>`
      + (SW.colorCol ? ` (grouped by <b>${SW.colorCol}</b>)` : "");
  body.innerHTML = `
    <div class="wiz-card">
      <div class="wiz-card-title">What we detected</div>
      <div class="wiz-desc">${desc}</div>
      ${SW.chartType !== "scatter" ? `
      <div class="wiz-groups">
        ${SW.groups.map(g=>`<span class="wiz-group-chip">${esc(g.label)} <em>n=${g.values.length}</em></span>`).join("")}
      </div>` : ""}
    </div>`;
}

// Step 1 — Normality
function _stepNormality(body) {
  if (SW.chartType === "scatter") {
    // For scatter just show distribution info
    const sx = shapiroWilk(SW.xVals), sy = shapiroWilk(SW.yVals);
    body.innerHTML = `
      <div class="wiz-card">
        <div class="wiz-card-title">Distribution check (Shapiro-Wilk)</div>
        ${_normRow(SW.xCol, SW.xVals.length, sx)}
        ${_normRow(SW.yCol, SW.yVals.length, sy)}
      </div>`;
    SW.normality = [
      { label:SW.xCol, n:SW.xVals.length, ...sx },
      { label:SW.yCol, n:SW.yVals.length, ...sy },
    ];
    return;
  }
  body.innerHTML = `
    <div class="wiz-card">
      <div class="wiz-card-title">Distribution check (Shapiro-Wilk per group)</div>
      ${SW.normality.map(r => _normRow(r.label, r.n, r)).join("")}
      <div class="wiz-note">${_normalitySummary()}</div>
    </div>`;
}

function _normRow(label, n, sw) {
  if (n < 3) return `<div class="norm-row"><span class="norm-lbl">${esc(label)}</span><span class="norm-badge grey">n &lt; 3 — skip</span></div>`;
  if (isNaN(sw?.W)) return `<div class="norm-row"><span class="norm-lbl">${esc(label)}</span><span class="norm-badge grey">—</span></div>`;
  const pass = sw.p > 0.05;
  return `<div class="norm-row">
    <span class="norm-lbl">${esc(label)} <em>(n=${n})</em></span>
    <span class="norm-detail">W = ${sw.W}, p = ${sw.p}</span>
    <span class="norm-badge ${pass?"green":"red"}">${pass?"Normal ✓":"Non-normal ✗"}</span>
  </div>`;
}

function _normalitySummary() {
  const allNormal = SW.normality.every(r=>r.p>0.05||r.n<3);
  if (allNormal) return "✓ All groups appear normally distributed — parametric tests are appropriate.";
  return "⚠ At least one group is non-normal — a non-parametric test is recommended.";
}

// Step 2 — Test selection
function _stepTestSelect(body) {
  const rec = _recommendTest();
  if (!SW.testName) SW.testName = rec;

  const opts = SW.chartType === "scatter"
    ? [
        { id:"pearson",    label:"Pearson correlation",  desc:"Linear relationship, normally distributed data" },
        { id:"spearman",   label:"Spearman correlation", desc:"Monotonic relationship, non-normal or ordinal data" },
      ]
    : SW.chartType === "twoGroup"
    ? [
        { id:"welch",       label:"Welch's t-test",     desc:"Two groups, continuous data, assumes normality" },
        { id:"mannwhitney", label:"Mann-Whitney U",      desc:"Two groups, non-parametric (no normality assumption)" },
      ]
    : [
        { id:"anova",    label:"One-way ANOVA",       desc:"3+ groups, assumes normality and equal variance" },
        { id:"kruskal",  label:"Kruskal-Wallis H",    desc:"3+ groups, non-parametric" },
      ];

  body.innerHTML = `
    <div class="wiz-card">
      <div class="wiz-card-title">Choose statistical test</div>
      <div class="wiz-test-opts">
        ${opts.map(o=>`
          <label class="wiz-test-opt ${o.id===rec?"wiz-recommended":""}">
            <input type="radio" name="wizTest" value="${o.id}" ${SW.testName===o.id?"checked":""}
                   onchange="SW.testName=this.value">
            <div>
              <div class="wiz-test-name">${o.label}${o.id===rec?' <span class="rec-badge">Recommended</span>':""}</div>
              <div class="wiz-test-desc">${o.desc}</div>
            </div>
          </label>`).join("")}
      </div>
    </div>`;
}

// Step 3 — Options (paired, tails, post-hoc, comparison style)
function _stepOptions(body) {
  const isTwoGroup = SW.chartType === "twoGroup";
  body.innerHTML = `
    <div class="wiz-card">
      <div class="wiz-card-title">Test options</div>

      ${isTwoGroup ? `
      <div class="wiz-opt-row">
        <span class="wiz-opt-lbl">Samples are</span>
        <div class="pill-row">
          <div class="pill ${!SW.paired?"on":""}" onclick="SW.paired=false;_syncPairedPills(this)">Unpaired</div>
          <div class="pill ${SW.paired?"on":""}"  onclick="SW.paired=true;_syncPairedPills(this)">Paired</div>
        </div>
      </div>
      <div class="wiz-opt-row">
        <span class="wiz-opt-lbl">Tails</span>
        <div class="pill-row">
          <div class="pill ${SW.tails===2?"on":""}" onclick="SW.tails=2;_syncTailPills(this)">Two-tailed</div>
          <div class="pill ${SW.tails===1?"on":""}" onclick="SW.tails=1;_syncTailPills(this)">One-tailed</div>
        </div>
      </div>` : `
      <div class="wiz-opt-row">
        <span class="wiz-opt-lbl">Post-hoc correction</span>
        <div class="pill-row" style="flex-wrap:wrap">
          <div class="pill ${SW.posthoc==="bonferroni"?"on":""}" onclick="SW.posthoc='bonferroni';_syncPosthocPills(this)">Bonferroni</div>
          <div class="pill ${SW.posthoc==="tukey"?"on":""}"      onclick="SW.posthoc='tukey';_syncPosthocPills(this)">Tukey HSD</div>
          <div class="pill ${SW.posthoc==="bh"?"on":""}"         onclick="SW.posthoc='bh';_syncPosthocPills(this)">BH / FDR</div>
          <div class="pill ${SW.posthoc==="none"?"on":""}"        onclick="SW.posthoc='none';_syncPosthocPills(this)">None</div>
        </div>
      </div>
      <div class="wiz-opt-row">
        <span class="wiz-opt-lbl">Compare</span>
        <div class="pill-row">
          <div class="pill ${SW.compStyle==="pairwise"?"on":""}"  onclick="SW.compStyle='pairwise';_syncCompPills(this);document.getElementById('refGroupRow').style.display='none'">All pairwise</div>
          <div class="pill ${SW.compStyle==="reference"?"on":""}" onclick="SW.compStyle='reference';_syncCompPills(this);document.getElementById('refGroupRow').style.display=''">vs. reference</div>
        </div>
      </div>
      <div id="refGroupRow" style="${SW.compStyle==="reference"?"":"display:none"};margin-top:6px">
        <span class="wiz-opt-lbl">Reference group</span>
        <select onchange="SW.refGroup=this.value" style="margin-top:4px;width:100%;padding:6px 8px;border:1px solid #dde1e7;border-radius:7px;font-size:13px">
          ${SW.groups.map(g=>`<option value="${esc(g.label)}" ${g.label===SW.refGroup?"selected":""}>${esc(g.label)}</option>`).join("")}
        </select>
      </div>`}
    </div>`;
}

function _syncPairedPills(el)  { el.closest(".pill-row").querySelectorAll(".pill").forEach(p=>p.classList.remove("on")); el.classList.add("on"); }
function _syncPosthocPills(el) { el.closest(".pill-row").querySelectorAll(".pill").forEach(p=>p.classList.remove("on")); el.classList.add("on"); }
function _syncTailPills(el)    { el.closest(".pill-row").querySelectorAll(".pill").forEach(p=>p.classList.remove("on")); el.classList.add("on"); }
function _syncCompPills(el)    { el.closest(".pill-row").querySelectorAll(".pill").forEach(p=>p.classList.remove("on")); el.classList.add("on"); }

// Step 4 — Confirm summary
function _stepConfirm(body) {
  const testLabels = {
    welch:"Welch's t-test", paired:"Paired t-test", mannwhitney:"Mann-Whitney U",
    anova:"One-way ANOVA", kruskal:"Kruskal-Wallis H",
    pearson:"Pearson correlation", spearman:"Spearman correlation",
  };
  const phLabels = { bonferroni:"Bonferroni", tukey:"Tukey HSD", bh:"BH/FDR", none:"None" };

  const rows = [
    ["Chart",        SW.chartType === "scatter" ? "Scatter (correlation)" : `${SW.groups.length} groups`],
    ["Test",         testLabels[SW.testName] || SW.testName],
    SW.chartType !== "scatter" && SW.chartType !== "twoGroup"
      ? ["Post-hoc",   phLabels[SW.posthoc]]  : null,
    SW.chartType !== "scatter" && SW.chartType !== "twoGroup"
      ? ["Compare",    SW.compStyle === "pairwise" ? "All pairwise" : `vs. ${SW.refGroup}`] : null,
    SW.chartType === "twoGroup"
      ? ["Paired",     SW.paired ? "Yes" : "No"] : null,
    SW.chartType === "twoGroup"
      ? ["Tails",      SW.tails + "-tailed"]    : null,
  ].filter(Boolean);

  body.innerHTML = `
    <div class="wiz-card">
      <div class="wiz-card-title">Confirm and run</div>
      <table class="wiz-summary-tbl">
        ${rows.map(([k,v])=>`<tr><td class="wiz-sum-k">${k}</td><td class="wiz-sum-v">${v}</td></tr>`).join("")}
      </table>
      <div class="wiz-note" style="margin-top:10px">Click <b>Run ▶</b> to compute and draw significance annotations.</div>
    </div>`;
}

// ── Run statistics ───────────────────────────────────────────────────
function _runStats() {
  let results = null;

  if (SW.chartType === "scatter") {
    const corr = SW.testName === "pearson"
      ? pearsonCorr(SW.xVals, SW.yVals)
      : spearmanCorr(SW.xVals, SW.yVals);
    const reg  = linearRegression(SW.xVals, SW.yVals);
    results = { type:"correlation", corr, reg, xCol:SW.xCol, yCol:SW.yCol };

  } else if (SW.chartType === "twoGroup") {
    const a = SW.groups[0].values, b = SW.groups[1].values;
    let testRes;
    if (SW.paired && SW.testName === "welch") SW.testName = "paired";
    if (SW.testName === "paired")      testRes = tTestPaired(a, b, SW.tails);
    else if (SW.testName === "welch")  testRes = tTestWelch(a, b, SW.tails);
    else                               testRes = mannWhitneyU(a, b, SW.tails);
    results = {
      type: "twoGroup",
      test: SW.testName, testRes,
      groups: SW.groups,
      desc: SW.groups.map(g=>({ label:g.label, ...descStats(g.values) })),
    };

  } else {
    // multiGroup
    let mainRes;
    if (SW.testName === "anova") mainRes = oneWayAnova(SW.groups);
    else                         mainRes = kruskalWallis(SW.groups);

    // Build pairs based on comparison style
    let pairs = SW.compStyle === "pairwise"
      ? _allPairs()
      : _refPairs();

    // Run pairwise tests
    pairs = _runPairwiseTests(pairs);

    // Apply post-hoc correction
    if (SW.posthoc !== "none" && SW.posthoc !== "tukey") {
      pairs = applyCorrection(pairs, SW.posthoc);
    } else if (SW.posthoc === "tukey" && SW.testName === "anova") {
      const tukeyPairs = postHocTukey(SW.groups);
      pairs = pairs.map(p => {
        const tp = tukeyPairs.find(t=>
          (t.groupA===p.groupA&&t.groupB===p.groupB)||(t.groupA===p.groupB&&t.groupB===p.groupA));
        return tp ? {...p, p_adj:tp.p_adj} : p;
      });
    } else if (SW.posthoc === "tukey" && SW.testName === "kruskal") {
      // Dunn's post-hoc
      const dunnPairs = postHocDunn(SW.groups);
      pairs = pairs.map(p => {
        const dp = dunnPairs.find(d=>
          (d.groupA===p.groupA&&d.groupB===p.groupB)||(d.groupA===p.groupB&&d.groupB===p.groupA));
        return dp ? {...p, p_raw:dp.p_raw, p_adj:dp.p_raw} : p;
      });
    }

    results = {
      type: "multiGroup",
      test: SW.testName, mainRes, pairs,
      groups: SW.groups,
      desc: SW.groups.map(g=>({ label:g.label, ...descStats(g.values) })),
    };
  }

  SW.results  = results;
  SW.applied  = true;
  closeStatsWizard();

  renderBrackets();
  _populateStatsTab();
  _generateStatsLegend();
}

function _allPairs() {
  const pairs = [];
  for(let i=0;i<SW.groups.length;i++)
    for(let j=i+1;j<SW.groups.length;j++)
      pairs.push({ groupA:SW.groups[i].label, groupB:SW.groups[j].label });
  return pairs;
}

function _refPairs() {
  return SW.groups
    .filter(g=>g.label!==SW.refGroup)
    .map(g=>({ groupA:SW.refGroup, groupB:g.label }));
}

function _runPairwiseTests(pairs) {
  return pairs.map(p => {
    const a = SW.groups.find(g=>g.label===p.groupA)?.values || [];
    const b = SW.groups.find(g=>g.label===p.groupB)?.values || [];
    let res;
    if (SW.testName === "anova") res = tTestWelch(a, b);
    else                         res = mannWhitneyU(a, b);
    const p_raw = res?.p ?? 1;
    return { ...p, p_raw:+p_raw.toFixed(6), p_adj:+p_raw.toFixed(6) };
  });
}

// ── Bracket rendering (Prism-style, paper x coords) ──────────────────
function renderBrackets() {
  if (!SW.results || !SW.applied) return;
  const res = SW.results;

  if (res.type === "correlation") {
    _renderRegressionLine(res);
    return;
  }

  // If a facet column is active, route to per-panel renderer
  const facetCol = document.getElementById("selFacet")?.value || null;
  if (facetCol) {
    _renderBracketsForFacets(facetCol);
    return;
  }

  const categories = SW.groups.map(g => g.label);
  const n = categories.length;
  if (n < 2) return;

  const pairs = res.type === "twoGroup"
    ? [{ groupA: SW.groups[0].label, groupB: SW.groups[1].label,
         p_adj: res.testRes?.p ?? 1, p_raw: res.testRes?.p ?? 1 }]
    : (res.pairs || []);

  // Y range from current layout
  const plotDiv = document.getElementById("plot");
  const layout  = plotDiv?._fullLayout;
  const yRange  = layout?.yaxis?.range;
  const allVals = SW.groups.flatMap(g => g.values).filter(v => isFinite(v));
  if (!allVals.length) return;
  const yMin  = yRange ? yRange[0] : Math.min(...allVals);
  const yMax  = yRange ? yRange[1] : Math.max(...allVals) * 1.1;
  const ySpan = Math.abs(yMax - yMin);

  // Map category index → paper x coordinate.
  // Plotly categorical axis: range is [-0.5, n-0.5]; categories evenly fill
  // the xaxis domain (defaults [0,1] for single-subplot charts).
  const xDomain = layout?.xaxis?.domain || [0, 1];
  const domW    = xDomain[1] - xDomain[0];
  const paperX  = i => xDomain[0] + ((i + 0.5) / n) * domW;

  // Sort: shorter span → lower bracket (wider spans rise higher, Prism style)
  const sorted = [...pairs].sort((a, b) => {
    const sa = Math.abs(categories.indexOf(a.groupA) - categories.indexOf(a.groupB));
    const sb = Math.abs(categories.indexOf(b.groupA) - categories.indexOf(b.groupB));
    return sa - sb;
  });

  const TICK_H = ySpan * 0.035;  // how far the vertical ticks drop below the bar
  const STEP   = ySpan * 0.09;   // vertical gap between successive brackets
  const BASE   = yMax + ySpan * 0.06;  // first bracket starts just above data max

  const shapes = [], annotations = [];

  sorted.forEach((p, idx) => {
    const iA = categories.indexOf(p.groupA);
    const iB = categories.indexOf(p.groupB);
    if (iA < 0 || iB < 0) return;

    const i0   = Math.min(iA, iB), i1 = Math.max(iA, iB);
    const x0p  = paperX(i0), x1p = paperX(i1);
    const xMid = (x0p + x1p) / 2;
    const yBar  = BASE + idx * STEP;
    const yTick = yBar - TICK_H;

    const useP = SW.posthoc === "none" ? p.p_raw : p.p_adj;
    const lbl  = fmtP(useP, SW.pDisplay);

    // Horizontal bar
    shapes.push({ type:"line", xref:"paper", yref:"y",
      x0:x0p, x1:x1p, y0:yBar, y1:yBar,
      line:{ color:"#333", width:1.5 } });
    // Left tick — drops DOWN from bar
    shapes.push({ type:"line", xref:"paper", yref:"y",
      x0:x0p, x1:x0p, y0:yTick, y1:yBar,
      line:{ color:"#333", width:1.5 } });
    // Right tick — drops DOWN from bar
    shapes.push({ type:"line", xref:"paper", yref:"y",
      x0:x1p, x1:x1p, y0:yTick, y1:yBar,
      line:{ color:"#333", width:1.5 } });
    // P-value label sits directly above the bar
    annotations.push({
      xref:"paper", yref:"y",
      x: xMid, y: yBar + ySpan * 0.012,
      text: lbl, showarrow: false,
      font:{ size: Math.max(9, ST.fontSize - 1), color:"#222" },
      xanchor:"center", yanchor:"bottom",
    });
  });

  const newYMax = BASE + sorted.length * STEP + STEP * 0.5;
  Plotly.relayout("plot", {
    shapes, annotations,
    "yaxis.range": [yMin - ySpan * 0.02, newYMax],
  });
}

function _clearBrackets() {
  // Guard: _fullLayout is only set after Plotly.react/newPlot has run at least once
  const plotDiv = document.getElementById("plot");
  if (!plotDiv || !plotDiv._fullLayout) return;
  try {
    // Preserve paper/paper annotations (facet panel titles); only remove bracket annotations
    const keepAnnots = (plotDiv._fullLayout.annotations || [])
      .filter(a => a.yref === "paper");
    Plotly.relayout("plot", { shapes:[], annotations: keepAnnots }).catch(()=>{});
  } catch(e) {}
}

// ── Per-facet bracket rendering ───────────────────────────────────────
// Called instead of renderBrackets() when a facet column is active.
// Re-runs the statistical test independently for each facet panel and
// draws Prism-style brackets using that panel's paper-coordinate domain.
function _renderBracketsForFacets(facetCol) {
  const plotDiv = document.getElementById("plot");
  if (!plotDiv?._fullLayout) return;

  const res = SW.results;
  const ct  = ST.chartType;

  // Same effective data used by the chart renderer
  const effectiveData = (FORMAT === "wide" && ct !== "heatmap") ? wideToLong() : DATA;

  // Facet values in same sorted order as renderFaceted()
  const facetVals = [...new Set(effectiveData.map(d => String(d[facetCol] ?? "(empty)")))]
    .sort((a, b) => a.localeCompare(b));

  // Preserve existing paper-space annotations (facet panel title labels)
  const existingAnnots = (plotDiv._fullLayout.annotations || [])
    .filter(a => a.yref === "paper");

  const allShapes    = [];
  const bracketAnnots = [];
  const relayoutUpd  = {};

  facetVals.forEach((fv, fi) => {
    const axId = fi === 0 ? "" : String(fi + 1);
    const xKey = `xaxis${axId}`;
    const yKey = `yaxis${axId}`;
    const yRef = `y${axId}`;

    const xDomain = plotDiv._fullLayout[xKey]?.domain || [0, 1];
    const yRange  = plotDiv._fullLayout[yKey]?.range;

    // --- Filter data to this facet panel ---
    const subset = effectiveData.filter(d => String(d[facetCol] ?? "(empty)") === fv);

    // Rebuild groups using same grouping key as the wizard
    const groupKey = SW.colorCol || SW.xCol;
    const labels   = [...new Set(subset.map(d => String(d[groupKey])))];
    const facetGroups = labels
      .map(lbl => ({
        label:  lbl,
        values: subset.filter(d => String(d[groupKey]) === lbl)
                      .map(d => parseFloat(d[SW.yCol])).filter(v => isFinite(v)),
      }))
      .filter(g => g.values.length > 0);

    if (facetGroups.length < 2) return;

    // --- Run the same test that the wizard chose, on this facet's groups ---
    let facetPairs;

    if (res.type === "twoGroup") {
      const a  = facetGroups[0]?.values || [];
      const b  = facetGroups[1]?.values || [];
      let tr;
      if      (SW.testName === "paired")      tr = tTestPaired(a, b, SW.tails);
      else if (SW.testName === "mannwhitney") tr = mannWhitneyU(a, b, SW.tails);
      else                                    tr = tTestWelch(a, b, SW.tails);
      const p = +(tr?.p ?? 1).toFixed(6);
      facetPairs = [{ groupA: facetGroups[0].label, groupB: facetGroups[1].label,
                      p_raw: p, p_adj: p }];

    } else {
      // multiGroup: build raw pairs respecting compStyle/refGroup
      const hasRef  = facetGroups.some(g => g.label === SW.refGroup);
      const rawPairs = (SW.compStyle === "reference" && hasRef)
        ? facetGroups.filter(g => g.label !== SW.refGroup)
                     .map(g => ({ groupA: SW.refGroup, groupB: g.label }))
        : (() => {
            const ps = [];
            for (let i = 0; i < facetGroups.length; i++)
              for (let j = i + 1; j < facetGroups.length; j++)
                ps.push({ groupA: facetGroups[i].label, groupB: facetGroups[j].label });
            return ps;
          })();

      // Run pairwise tests (use appropriate test per wizard choice)
      facetPairs = rawPairs.map(p => {
        const a = facetGroups.find(g => g.label === p.groupA)?.values || [];
        const b = facetGroups.find(g => g.label === p.groupB)?.values || [];
        let tr;
        if (SW.testName === "kruskal") tr = mannWhitneyU(a, b);
        else                           tr = tTestWelch(a, b);
        const p_raw = +(tr?.p ?? 1).toFixed(6);
        return { ...p, p_raw, p_adj: p_raw };
      });

      // Apply multiple-testing correction
      if (SW.posthoc !== "none" && facetPairs.length > 1) {
        facetPairs = applyCorrection(facetPairs, SW.posthoc);
      }
    }

    if (!facetPairs.length) return;

    // --- Geometry: Y space and paper X mapping ---
    const allVals = facetGroups.flatMap(g => g.values).filter(v => isFinite(v));
    if (!allVals.length) return;
    const yMin  = yRange ? yRange[0] : Math.min(...allVals);
    const yMax  = yRange ? yRange[1] : Math.max(...allVals) * 1.1;
    const ySpan = Math.abs(yMax - yMin);

    const categories = facetGroups.map(g => g.label);
    const n     = categories.length;
    const domW  = xDomain[1] - xDomain[0];
    const paperX = i => xDomain[0] + ((i + 0.5) / n) * domW;

    // Sort: shorter span → lower bracket (Prism-style stacking)
    const sorted = [...facetPairs].sort((a, b) => {
      const sa = Math.abs(categories.indexOf(a.groupA) - categories.indexOf(a.groupB));
      const sb = Math.abs(categories.indexOf(b.groupA) - categories.indexOf(b.groupB));
      return sa - sb;
    });

    const TICK_H = ySpan * 0.035;
    const STEP   = ySpan * 0.09;
    const BASE   = yMax + ySpan * 0.06;

    sorted.forEach((p, idx) => {
      const iA = categories.indexOf(p.groupA);
      const iB = categories.indexOf(p.groupB);
      if (iA < 0 || iB < 0) return;

      const i0   = Math.min(iA, iB), i1 = Math.max(iA, iB);
      const x0p  = paperX(i0), x1p = paperX(i1);
      const xMid = (x0p + x1p) / 2;
      const yBar  = BASE + idx * STEP;
      const yTick = yBar - TICK_H;

      const useP = SW.posthoc === "none" ? p.p_raw : p.p_adj;
      const lbl  = fmtP(useP, SW.pDisplay);

      // Horizontal bar
      allShapes.push({ type:"line", xref:"paper", yref:yRef,
        x0:x0p, x1:x1p, y0:yBar, y1:yBar, line:{color:"#333",width:1.5} });
      // Left tick (drops DOWN from bar)
      allShapes.push({ type:"line", xref:"paper", yref:yRef,
        x0:x0p, x1:x0p, y0:yTick, y1:yBar, line:{color:"#333",width:1.5} });
      // Right tick (drops DOWN from bar)
      allShapes.push({ type:"line", xref:"paper", yref:yRef,
        x0:x1p, x1:x1p, y0:yTick, y1:yBar, line:{color:"#333",width:1.5} });
      // P-value label above bar
      bracketAnnots.push({
        xref:"paper", yref:yRef,
        x:xMid, y:yBar + ySpan * 0.012,
        text:lbl, showarrow:false,
        font:{ size:Math.max(9, ST.fontSize-1), color:"#222" },
        xanchor:"center", yanchor:"bottom",
      });
    });

    // Expand this panel's Y axis to accommodate all brackets
    const newYMax = BASE + sorted.length * STEP + STEP * 0.5;
    relayoutUpd[`${yKey}.range`] = [yMin - ySpan * 0.02, newYMax];
  });

  // Merge bracket annotations with preserved facet title annotations
  relayoutUpd.shapes      = allShapes;
  relayoutUpd.annotations = [...existingAnnots, ...bracketAnnots];

  Plotly.relayout("plot", relayoutUpd);
}

function _renderRegressionLine(res) {
  if (!res.reg) return;
  const { slope, intercept, r2 } = res.reg;
  const xSorted = [...SW.xVals].sort((a,b)=>a-b);
  const xLine   = [xSorted[0], xSorted[xSorted.length-1]];
  const yLine   = xLine.map(x=>slope*x+intercept);

  const corrLabel = res.corr
    ? `${res.corr.method}: r = ${res.corr.r}, ${fmtP(res.corr.p, SW.pDisplay)}`
    : "";
  const regLabel = `y = ${slope.toFixed(3)}x + ${intercept.toFixed(3)}, R² = ${r2}`;

  Plotly.addTraces("plot", [{
    type:"scatter", mode:"lines",
    x: xLine, y: yLine,
    name: "Fit",
    line:{ color:"#e15759", width:2, dash:"dash" },
    hoverinfo:"skip",
  }]);

  Plotly.relayout("plot", {
    annotations:[{
      xref:"paper", yref:"paper", x:0.02, y:0.97,
      text: corrLabel + "<br>" + regLabel,
      showarrow:false, xanchor:"left", yanchor:"top",
      font:{ size:Math.max(10,ST.fontSize-1), color:"#333" },
      bgcolor:"rgba(255,255,255,0.85)",
      bordercolor:"#ccc", borderwidth:1,
    }],
  });
}

// ── Stats tab population ──────────────────────────────────────────────
function _populateStatsTab() {
  const panel = document.getElementById("statsTabContent");
  if (!panel || !SW.results) return;

  const res = SW.results;
  let html = "";

  // ── Main test result ──
  if (res.type === "correlation") {
    const c = res.corr, r = res.reg;
    html += `<div class="stats-section-head">Correlation</div>
      <table class="stats-tbl">
        <tr><td>Method</td><td>${c?.method || "—"}</td></tr>
        <tr><td>r</td><td>${c?.r ?? "—"}</td></tr>
        <tr><td>R²</td><td>${c?.r2 ?? "—"}</td></tr>
        <tr><td>p-value</td><td>${c ? fmtP(c.p,"exact") : "—"}</td></tr>
        <tr><td>df</td><td>${c?.df ?? "—"}</td></tr>
        ${r ? `<tr><td>Slope</td><td>${r.slope}</td></tr>
               <tr><td>Intercept</td><td>${r.intercept}</td></tr>` : ""}
      </table>`;

  } else if (res.type === "twoGroup") {
    const t = res.testRes;
    const testLbl = {welch:"Welch's t-test",paired:"Paired t-test",mannwhitney:"Mann-Whitney U"}[res.test]||res.test;
    html += `<div class="stats-section-head">${testLbl}</div>
      <table class="stats-tbl">
        ${"t" in (t||{}) ? `<tr><td>t</td><td>${t.t}</td></tr><tr><td>df</td><td>${t.df}</td></tr>` : ""}
        ${"U" in (t||{}) ? `<tr><td>U</td><td>${t.U}</td></tr><tr><td>z</td><td>${t.z}</td></tr>` : ""}
        <tr><td>p-value</td><td>${t ? fmtP(t.p,"exact") : "—"}</td></tr>
        ${"cohensD" in (t||{}) ? `<tr><td>Cohen's d</td><td>${t.cohensD}</td></tr>` : ""}
        ${"r" in (t||{}) ? `<tr><td>Effect r</td><td>${t.r}</td></tr>` : ""}
      </table>`;

  } else {
    const m = res.mainRes;
    const testLbl = {anova:"One-way ANOVA",kruskal:"Kruskal-Wallis H"}[res.test]||res.test;
    html += `<div class="stats-section-head">${testLbl}</div>
      <table class="stats-tbl">
        ${"F" in (m||{}) ? `<tr><td>F (${m.dfb}, ${m.dfw})</td><td>${m.F}</td></tr>` : ""}
        ${"H" in (m||{}) ? `<tr><td>H (df=${m.df})</td><td>${m.H}</td></tr>` : ""}
        <tr><td>p-value</td><td>${m ? fmtP(m.p,"exact") : "—"}</td></tr>
        ${"eta2" in (m||{}) ? `<tr><td>η²</td><td>${m.eta2}</td></tr>` : ""}
      </table>`;

    // Pairwise table
    if (res.pairs?.length) {
      const phLbl = {bonferroni:"Bonferroni",tukey:"Tukey HSD",bh:"BH/FDR",none:"No correction"}[SW.posthoc]||SW.posthoc;
      html += `<div class="stats-section-head" style="margin-top:12px">Post-hoc comparisons (${phLbl})</div>
        <table class="stats-tbl stats-pairs-tbl">
          <thead><tr><th>Group A</th><th>Group B</th><th>p (raw)</th><th>p (adj)</th><th></th></tr></thead>
          <tbody>
            ${res.pairs.map(p=>`
              <tr>
                <td>${esc(p.groupA)}</td><td>${esc(p.groupB)}</td>
                <td>${fmtP(p.p_raw,"exact")}</td>
                <td>${fmtP(p.p_adj,"exact")}</td>
                <td class="sig-stars">${sigStars(p.p_adj)}</td>
              </tr>`).join("")}
          </tbody>
        </table>`;
    }
  }

  // ── Descriptive stats table ──
  if (res.desc?.length) {
    html += `<div class="stats-section-head" style="margin-top:12px">Descriptive statistics</div>
      <div style="overflow-x:auto">
      <table class="stats-tbl stats-desc-tbl">
        <thead><tr><th>Group</th><th>N</th><th>Mean</th><th>SD</th><th>SEM</th><th>Median</th><th>Q1</th><th>Q3</th><th>Min</th><th>Max</th></tr></thead>
        <tbody>
          ${res.desc.map(d=>`
            <tr>
              <td>${esc(d.label)}</td>
              <td>${d.n}</td>
              <td>${_fmt(d.mean)}</td><td>${_fmt(d.sd)}</td><td>${_fmt(d.sem)}</td>
              <td>${_fmt(d.median)}</td><td>${_fmt(d.q1)}</td><td>${_fmt(d.q3)}</td>
              <td>${_fmt(d.min)}</td><td>${_fmt(d.max)}</td>
            </tr>`).join("")}
        </tbody>
      </table></div>`;
  }

  // Export button
  html += `<button class="btn btn-grey" style="margin-top:12px;width:100%" onclick="exportStatsCSV()">⬇ Export stats CSV</button>`;

  panel.innerHTML = html;
}

function _fmt(v) {
  if (v===undefined||v===null||isNaN(v)) return "—";
  return Math.abs(v)<0.001||Math.abs(v)>=10000 ? v.toExponential(3) : +v.toPrecision(4)+"";
}

function _clearStatsTab() {
  const panel = document.getElementById("statsTabContent");
  if (panel) panel.innerHTML = `<div class="stats-empty">Run statistics to see results here.</div>`;
}

// Export stats as CSV
function exportStatsCSV() {
  if (!SW.results?.desc) return;
  const res = SW.results;
  const rows = [["Group","N","Mean","SD","SEM","Median","Q1","Q3","Min","Max"]];
  res.desc.forEach(d=>{
    rows.push([d.label,d.n,d.mean,d.sd,d.sem,d.median,d.q1,d.q3,d.min,d.max]);
  });
  if (res.pairs) {
    rows.push([]);
    rows.push(["Group A","Group B","p (raw)","p (adj)","Significance"]);
    res.pairs.forEach(p=>rows.push([p.groupA,p.groupB,p.p_raw,p.p_adj,sigStars(p.p_adj)]));
  }
  const csv = rows.map(r=>r.map(c=>`"${String(c??"")}"`).join(",")).join("\n");
  const a = document.createElement("a");
  a.href = "data:text/csv;charset=utf-8,"+encodeURIComponent(csv);
  a.download = "stats_results.csv";
  a.click();
}

// ── Figure legend — Section 1 (always generated after render) ────────
function _generateFigureLegend() {
  const ta = document.getElementById("legendFigureTextarea");
  if (!ta) return;

  const ct = ST.chartType;
  const xCol    = document.getElementById("selX")?.value    || "[X]";
  const yCol    = document.getElementById("selY")?.value    || "[Y]";
  const colorCol= document.getElementById("selColor")?.value || "";

  const chartDesc = {
    bar:"Bar graph", scatter:"Scatter plot", line:"Line graph",
    box:"Box plot", violin:"Violin plot", histogram:"Histogram",
    heatmap:"Heatmap", alluvial:"Alluvial diagram",
  }[ct] || "Figure";

  let desc = `Figure X. ${chartDesc}`;

  if (ct === "heatmap") {
    if (FORMAT === "wide") {
      const rowSplit  = document.getElementById("selRowSplit")?.value || "";
      const rowLabel  = document.getElementById("selRowLabel")?.value || "features";
      const transforms = [];
      if (ST.log2)   transforms.push("log₂(x+1) transformed");
      if (ST.center) transforms.push("row-mean centered");
      desc += ` showing ${rowLabel} expression across [samples]`;
      if (rowSplit) desc += `, annotated by ${rowSplit}`;
      if (transforms.length) desc += `. Values are ${transforms.join(" and ")}`;
      desc += `. Color scale: ${ST.center ? "diverging (centered on zero)" : "sequential (low→high)"}.`;
    } else {
      desc += ` of ${colorCol || "[Z]"} values for ${yCol} (rows) × ${xCol} (columns).`;
    }
  } else if (ct === "histogram") {
    desc += ` of ${xCol}` + (colorCol ? `, colored by ${colorCol}` : "") + ".";
  } else if (ct === "scatter") {
    desc += ` of ${yCol} vs. ${xCol}` + (colorCol ? `, colored by ${colorCol}` : "") + ".";
  } else if (ct === "alluvial") {
    desc += ` showing flow from ${xCol}` + (colorCol ? ` through ${colorCol}` : "") + ` to ${yCol}.`;
  } else {
    // bar, box, violin, line
    desc += ` of ${yCol} by ${xCol}` + (colorCol ? `, grouped by ${colorCol}` : "") + ".";
    const statDesc = {
      bar:    "Bars show mean ± SEM.",
      box:    "Box plots show median with interquartile range; whiskers extend to 1.5× IQR.",
      violin: "Violin plots show kernel density estimate with embedded box (median, IQR) and mean line.",
      line:   "Lines connect mean values.",
    }[ct];
    if (statDesc) desc += " " + statDesc;

    // Group N summary — wide and long format handled separately
    const effectiveData = (FORMAT === "wide" && ct !== "heatmap") ? wideToLong() : DATA;
    const facetActive   = !!(document.getElementById("selFacet")?.value);
    const nSuffix       = facetActive ? " (per facet)" : "";
    if (effectiveData.length) {
      if (FORMAT === "wide") {
        // For wide format n always counts unique samples (not gene×sample rows).
        // Group by colorCol only (xCol can be a high-cardinality gene list).
        if (colorCol) {
          const colorGroups = [...new Set(effectiveData.map(d => String(d[colorCol])))];
          if (colorGroups.length <= 10) {
            const summary = colorGroups.map(g => {
              const gr = effectiveData.filter(d => String(d[colorCol]) === g);
              const n  = new Set(gr.map(d => d.Sample)).size;
              return `${g} (n=${n} samples${nSuffix})`;
            });
            desc += " " + summary.join(", ") + ".";
          } else {
            const n = new Set(effectiveData.map(d => d.Sample)).size;
            desc += ` n=${n} samples total.`;
          }
        } else {
          const n = new Set(effectiveData.map(d => d.Sample)).size;
          desc += ` n=${n} samples${nSuffix}.`;
        }
      } else {
        // Long format: enumerate groups with observation counts
        const groupKey = colorCol || xCol;
        if (groupKey) {
          const groups = [...new Set(effectiveData.map(d => String(d[groupKey])))];
          if (groups.length <= 10) {
            const summary = groups.map(g => {
              const gr = effectiveData.filter(d => String(d[groupKey]) === g);
              const n  = gr.map(d => parseFloat(d[yCol])).filter(v => isFinite(v)).length;
              return `${g} (n=${n}${nSuffix})`;
            });
            desc += " " + summary.join(", ") + ".";
          } else {
            desc += ` N = ${effectiveData.length} total observations.`;
          }
        }
      }
    }
  }

  ta.value = desc;
  // Always show the legend panel once we have data
  document.getElementById("legendPanel").style.display = "";
}

// ── Stats legend — Section 2 (generated after stats run) ─────────────
function _generateStatsLegend() {
  const ta = document.getElementById("legendStatsTextarea");
  if (!ta || !SW.results) return;

  const res = SW.results;
  let text = "";

  if (res.type === "twoGroup") {
    const testLbl = {
      welch:       "Welch's unpaired t-test",
      paired:      "paired Student's t-test",
      mannwhitney: "Mann-Whitney U test",
    }[res.test] || res.test;
    const tails = SW.tails === 1 ? "one-tailed" : "two-tailed";
    text = `Statistical comparison by ${tails} ${testLbl}`;
    if (res.testRes?.cohensD !== undefined) text += ` (Cohen's d = ${res.testRes.cohensD})`;
    text += ".";

  } else if (res.type === "multiGroup") {
    const testLbl = { anova:"one-way ANOVA", kruskal:"Kruskal-Wallis H test" }[res.test] || res.test;
    const phLbl   = { bonferroni:"Bonferroni", tukey:"Tukey HSD", bh:"Benjamini-Hochberg FDR", none:"no" }[SW.posthoc];
    const compDesc = SW.compStyle === "pairwise"
      ? "all pairwise comparisons"
      : `comparisons vs. ${SW.refGroup}`;
    text = `Statistical comparisons by ${testLbl} with post-hoc ${phLbl} correction for ${compDesc}`;
    if (res.mainRes?.eta2 !== undefined) text += ` (η² = ${res.mainRes.eta2})`;
    if (res.mainRes?.p    !== undefined) text += `; overall ${fmtP(res.mainRes.p,"exact")}`;
    text += ".";

  } else if (res.type === "correlation") {
    const method = res.corr?.method || "Pearson";
    text = `${method} correlation: r = ${res.corr?.r ?? "—"}, ${fmtP(res.corr?.p,"exact")}, R² = ${res.corr?.r2 ?? "—"}.`;
    if (res.reg) text += ` Regression line: y = ${res.reg.slope}x + ${res.reg.intercept}.`;
  }

  const sigKey = SW.pDisplay === "stars"
    ? " * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001; ns = not significant."
    : " Exact p-values shown above comparison brackets.";
  text += sigKey;

  ta.value = text;
  ta.removeAttribute("readonly");
}

function copyLegend() {
  const fig   = document.getElementById("legendFigureTextarea")?.value || "";
  const stats = document.getElementById("legendStatsTextarea")?.value  || "";
  const combined = [fig, stats].filter(Boolean).join(" ");
  navigator.clipboard?.writeText(combined).catch(() => {
    // Fallback for older browsers
    const ta = document.createElement("textarea");
    ta.value = combined;
    document.body.appendChild(ta);
    ta.select();
    document.execCommand("copy");
    document.body.removeChild(ta);
  });
  const btn = document.getElementById("copyLegendBtn");
  btn.textContent = "Copied ✓";
  setTimeout(() => { btn.textContent = "Copy All"; }, 2000);
}

function togglePDisplay(mode) {
  SW.pDisplay = mode;
  document.querySelectorAll("#pDisplayToggle .pill").forEach(p =>
    p.classList.toggle("on", p.dataset.pd === mode));
  if (SW.applied) {
    renderBrackets();
    _generateStatsLegend();
  }
}
