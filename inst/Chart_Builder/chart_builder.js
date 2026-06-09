/*! © 2024–2026 Nourhan Abdelfattah — scSidekick R package
 *  Bundled for use via scSidekick::ChartBuilder() only.
 *  Not licensed for standalone redistribution.
 */
// ═══════════════════════════════════════════════════════════════════
//  CHART RENDERERS + EXPORT  (depends on: config.js, statistics.js, ui_controls.js)
// ═══════════════════════════════════════════════════════════════════

function getColors()   {
  if (ST.catPal === "custom") return ST.customColors;
  return CAT_PALS[ST.catPal]?.colors || CAT_PALS.tableau.colors;
}
// Returns the color for a group — checks per-group overrides first, then palette
function getGroupColor(groupName, paletteIndex) {
  if (ST.groupColors && ST.groupColors[groupName]) return ST.groupColors[groupName];
  const colors = getColors();
  return colors[paletteIndex % colors.length];
}
// Color for jitter/dot overlay points
function pointColor(groupName, paletteIndex) {
  return ST.blackPoints ? "#333333" : getGroupColor(groupName, paletteIndex);
}
function hexToRgba(hex, alpha) {
  if (!hex || !hex.startsWith("#") || hex.length < 7) return `rgba(74,143,232,${alpha})`;
  const r = parseInt(hex.slice(1,3),16), g = parseInt(hex.slice(3,5),16), b = parseInt(hex.slice(5,7),16);
  return `rgba(${r},${g},${b},${alpha})`;
}
function getSeqScale() { return SEQ_PALS[ST.seqPal]?.scale  || "Viridis"; }
function getDivScale() { return DIV_PALS[ST.divPal]?.scale  || DIV_PALS.rbu.scale; }
// Active heatmap colorscale — follows the user's palette choice (seq or div),
// independent of the Center transform.
function getHeatScale() { return ST.heatScaleType === "div" ? getDivScale() : getSeqScale(); }
function gridColor()   { return ST.grid ? "#e8edf2" : "rgba(0,0,0,0)"; }
function yLabel(yCol, aggFn) {
  if (aggFn === "count") return "Count";
  if (aggFn === "none")  return yCol;
  return `${aggFn}(${yCol})`;
}

// ═══════════════════════════════════════════════════════════════════
//  MASTER RENDER DISPATCHER
// ═══════════════════════════════════════════════════════════════════
function render() {
  if (!DATA.length) return;
  updateAxisUI();

  const ct      = ST.chartType;
  const isWide  = FORMAT === "wide";
  const isWideHeat = isWide && ct === "heatmap";

  // Problem 5: for wide + non-heatmap, melt the data first
  let renderData = (isWide && !isWideHeat) ? wideToLong() : DATA;
  // Apply long-format column value filters
  if (!isWide && Object.keys(ST.colFilters || {}).length > 0) {
    renderData = renderData.filter(row =>
      Object.entries(ST.colFilters).every(([col, allowed]) =>
        !allowed || allowed.has(String(row[col] ?? ""))
      )
    );
  }

  let res;
  if (isWideHeat) {
    res = renderWideHeatmap();
  } else {
    // Alluvial can't be faceted — always ignore the facet selector for it
    const facetCol = (ct === "alluvial") ? null : (document.getElementById("selFacet")?.value || null);
    const nCols    = parseInt(document.getElementById("facetNCols")?.value) || ST.nFacetCols;
    const freeY    = ST.freeY;

    const rendererMap = {
      bar:       d => renderBar(d),
      scatter:   d => renderScatter(d),
      line:      d => renderLine(d),
      histogram: d => renderHistogram(d),
      box:       d => renderBox(d),
      violin:    d => renderViolin(d),
      heatmap:   d => renderHeatmap(d),
      alluvial:  d => renderAlluvial(d),
    };
    const baseRenderer = rendererMap[ct];
    if (!baseRenderer) return;

    if (facetCol) {
      res = renderFaceted(baseRenderer, renderData, facetCol, nCols, freeY);
    } else {
      res = baseRenderer(renderData);
    }
  }
  if (!res) return;

  const layout = {
    ...res.layout,
    paper_bgcolor: "white",
    plot_bgcolor:  ST.grid ? "#fafbfc" : "white",
    font: { family: FONT, size: ST.fontSize },
    margin: { l:70, r:40, t:30, b:80, ...(res.layout?.margin||{}) },
    legend: { ...(res.layout?.legend||{}), font:{size: Math.max(9, ST.fontSize - 1)} },
  };

  // Apply axis line styling (skip for alluvial/sankey which has no xaxis/yaxis)
  if (ct !== "alluvial") {
    const axisStyle = {
      linecolor:  ST.axisColor,
      linewidth:  ST.axisThick,
      showline:   ST.axisThick > 0,
    };
    const yRangeStyle = ST.yFromZero ? { rangemode: "tozero" } : {};
    // Ensure xaxis/yaxis exist, then apply to all xy-axis keys (covers facet sub-axes too)
    if (!layout.xaxis) layout.xaxis = {};
    if (!layout.yaxis) layout.yaxis = {};
    Object.keys(layout).forEach(k => {
      if (/^[xy]axis\d*$/.test(k)) {
        layout[k] = { ...layout[k], ...axisStyle };
        if (/^yaxis\d*$/.test(k)) layout[k] = { ...layout[k], ...yRangeStyle };
      }
    });
  }

  Plotly.react("plot", res.traces, layout, {
    responsive: true,
    displayModeBar: true,
    modeBarButtonsToRemove: ["select2d","lasso2d","autoScale2d","toggleSpikelines"],
    displaylogo: false,
    toImageButtonOptions: { format:"png", width:ST.exportW, height:ST.exportH, scale:1, filename:"chart" },
  }).then(() => {
    // Re-apply significance brackets after every re-render
    if (typeof SW !== "undefined" && SW.applied && SW.results) {
      if (SW.results.type !== "correlation") renderBrackets();
    }
    // Always refresh the figure description section of the legend
    try {
      if (typeof _generateFigureLegend === "function") _generateFigureLegend();
    } catch(e) { console.error("_generateFigureLegend error:", e); }
  });
}

// ═══════════════════════════════════════════════════════════════════
//  PROBLEM 1: GENERIC FACET WRAPPER
// ═══════════════════════════════════════════════════════════════════
function renderFaceted(baseRendererFn, data, facetCol, nCols, freeY) {
  const facetVals = [...new Set(data.map(d => String(d[facetCol] ?? "(empty)")))]
    .sort((a,b) => a.localeCompare(b));
  const nPanels = facetVals.length;
  if (!nPanels) return null;

  const nC = Math.min(Math.max(1, nCols), nPanels);
  const nR = Math.ceil(nPanels / nC);

  const GAP_X = 0.06, GAP_Y = 0.16;
  const panW  = (1 - (nC - 1) * GAP_X) / nC;
  const panH  = (1 - (nR - 1) * GAP_Y) / nR;

  // Compute global Y range for shared scale
  let globalYMin = Infinity, globalYMax = -Infinity;
  if (!freeY) {
    facetVals.forEach(fv => {
      const subset = data.filter(d => String(d[facetCol] ?? "(empty)") === fv);
      const res = baseRendererFn(subset);
      if (!res) return;
      res.traces.forEach(t => {
        (t.y || []).forEach(y => {
          if (y != null && !isNaN(+y)) {
            globalYMin = Math.min(globalYMin, +y);
            globalYMax = Math.max(globalYMax, +y);
          }
        });
      });
    });
    // Pad range
    const pad = (globalYMax - globalYMin) * 0.08 || 1;
    globalYMin -= pad; globalYMax += pad;
  }

  const allTraces = [], allAnnots = [], layoutAxes = {};

  facetVals.forEach((fv, fi) => {
    const col = fi % nC;
    const row = Math.floor(fi / nC);
    const x0  = col * (panW + GAP_X);
    const y1  = 1 - row * (panH + GAP_Y);
    const xd  = [x0, x0 + panW];
    const yd  = [y1 - panH, y1];

    const subset = data.filter(d => String(d[facetCol] ?? "(empty)") === fv);
    const res = baseRendererFn(subset);
    if (!res) return;

    // Axis ID strings: panel 0 → "x"/"y", panel N → "x(N+1)"/"y(N+1)"
    const axId  = fi === 0 ? "" : String(fi + 1);
    const xRef  = `x${axId}`;
    const yRef  = `y${axId}`;
    const xKey  = `xaxis${axId}`;
    const yKey  = `yaxis${axId}`;

    res.traces.forEach((t, ti) => {
      const remapped = { ...t };
      remapped.xaxis = xRef;
      remapped.yaxis = yRef;
      // Suppress colorbar on all but last panel
      if (remapped.colorbar && fi < nPanels - 1) remapped.colorbar = undefined;
      if (remapped.showscale !== undefined && fi < nPanels - 1) remapped.showscale = false;
      // Only show legend entries from the first panel
      if (fi > 0) remapped.showlegend = false;
      allTraces.push(remapped);
    });

    const srcX = res.layout?.xaxis || {};
    const srcY = res.layout?.yaxis || {};

    layoutAxes[xKey] = {
      ...srcX,
      domain: xd,
      anchor: yRef,
      showticklabels: row === nR - 1,
    };
    layoutAxes[yKey] = {
      ...srcY,
      domain: yd,
      anchor: xRef,
      showticklabels: col === 0,
      ...(!freeY && isFinite(globalYMin) ? { range: [globalYMin, globalYMax] } : {}),
    };

    allAnnots.push({
      text: `<b>${fv}</b>`,
      xref:"paper", yref:"paper",
      x: (xd[0] + xd[1]) / 2,
      y: yd[1] + 0.01,
      xanchor:"center", yanchor:"bottom",
      showarrow: false,
      font: { size: ST.fontSize, color:"#1e2a38", family:FONT },
    });
  });

  return {
    traces: allTraces,
    layout: {
      annotations: allAnnots,
      margin: { l:60, r:50, t:60, b:80 },
      ...layoutAxes,
    },
  };
}

// ── Bar ──────────────────────────────────────────────────────────────
function renderBar(data = DATA) {
  const xCol     = document.getElementById("selX").value;
  const yCol     = document.getElementById("selY").value;
  const colorCol = document.getElementById("selColor").value || null;
  const aggFn    = document.getElementById("selAgg").value;
  const sort     = document.getElementById("selSort").value;
  if (!xCol || !yCol) return;

  const agg    = aggregate(data, xCol, yCol, colorCol, aggFn);
  const colors = getColors();
  const xOrder = uniqueXOrder(agg, sort);
  const showMean = aggFn === "mean";

  // Helper: build error_y for one set of aggregate rows (ordered by xOrder)
  const errArr = (rows) => xOrder.map(xv => {
    const f = rows.find(d => d.x === xv);
    if (!f) return null;
    return ST.errorBars === "sd" ? f.sd : f.sem;
  });
  const errY = (rows) => showMean && ST.errorBars !== "none"
    ? { type:"data", array:errArr(rows), visible:true,
        color:"#444", thickness:1.5, width:5, symmetric:true }
    : { visible:false };

  let traces;
  if (colorCol) {
    const groups = [...new Set(agg.map(d => d.color))];

    // Detect diagonal pattern: each color group appears at exactly one x-position
    const colorXSets = new Map();
    agg.forEach(d => {
      if (!colorXSets.has(d.color)) colorXSets.set(d.color, new Set());
      colorXSets.get(d.color).add(d.x);
    });
    const isDiagonal = groups.every(g => (colorXSets.get(g)?.size || 0) <= 1);

    if (isDiagonal) {
      // One bar trace with per-bar color array, plus dummy scatter traces for the legend
      const colorArr = xOrder.map(xv => {
        const d = agg.find(a => a.x === xv);
        if (!d) return getColors()[0];
        const gi = groups.indexOf(d.color);
        return getGroupColor(String(d.color), gi >= 0 ? gi : 0);
      });
      const ys = xOrder.map(xv => { const f = agg.find(a => a.x === xv); return f ? f.y : null; });
      traces = [{
        type:"bar",
        x: xOrder, y: ys,
        error_y: errY(agg),
        marker:{ color:colorArr, opacity:.85, line:{color:"black",width:0.8} },
        text:      ST.labels ? ys.map(v => v!=null ? (+v.toFixed(2)).toString() : "") : [],
        textposition: ST.labels ? "outside" : "none",
        showlegend: false,
      }];
      // Dummy scatter traces for legend (one per color group)
      groups.forEach((g, gi) => {
        traces.push({
          type:"scatter", mode:"markers", name:String(g),
          x:[null], y:[null],
          showlegend:true,
          marker:{ color:getGroupColor(String(g), gi), size:11, symbol:"square" },
          legendrank: gi + 100,
        });
      });
      // Jitter point overlays for diagonal — attach to each x-group
      if (ST.showPoints !== "none") {
        groups.forEach((g, gi) => {
          const gd = data.filter(d => String(d[colorCol]) === g);
          const pts = gd.map(d => ({ x:String(d[xCol]), y:parseFloat(d[yCol]) }))
                        .filter(d => isFinite(d.y));
          if (!pts.length) return;
          traces.push({
            type:"box", name:String(g),
            x:pts.map(d=>d.x), y:pts.map(d=>d.y),
            boxpoints:"all", jitter:0.35, pointpos:0,
            whiskerwidth:0, boxmean:false,
            fillcolor:"rgba(0,0,0,0)", line:{color:"rgba(0,0,0,0)",width:0},
            marker:{ color:pointColor(String(g), gi), size:5, opacity:0.7,
                     line:{color:"white",width:0.5} },
            showlegend:false,
          });
        });
      }
    } else {
      // Normal grouped bar — use null for missing color×x combos
      traces = groups.map((g, gi) => {
        const gd = agg.filter(d => d.color === g);
        const ys = xOrder.map(xv => { const f=gd.find(d=>d.x===xv); return f ? f.y : null; });
        return {
          type:"bar", name:String(g),
          x:xOrder, y:ys,
          error_y: errY(gd),
          marker:{ color:getGroupColor(String(g), gi), opacity:.85, line:{color:"black",width:0.8} },
          text:      ST.labels ? ys.map(v => v!=null ? (+v.toFixed(2)).toString() : "") : [],
          textposition: ST.labels ? "outside" : "none",
        };
      });
      // Jitter point overlays (one invisible-box per color group)
      if (ST.showPoints !== "none") {
        groups.forEach((g, gi) => {
          const gd = data.filter(d => d[colorCol] === g);
          const pts = gd.map(d => ({ x:String(d[xCol]), y:parseFloat(d[yCol]) }))
                        .filter(d => isFinite(d.y));
          if (!pts.length) return;
          traces.push({
            type:"box", name:String(g),
            x:pts.map(d=>d.x), y:pts.map(d=>d.y),
            boxpoints:"all", jitter:0.35, pointpos:0,
            whiskerwidth:0, boxmean:false,
            fillcolor:"rgba(0,0,0,0)", line:{color:"rgba(0,0,0,0)",width:0},
            marker:{ color:pointColor(String(g), gi), size:5, opacity:0.7,
                     line:{color:"white",width:0.5} },
            showlegend:false,
          });
        });
      }
    }
  } else {
    const sorted = sortAgg(agg, sort);
    traces = [{
      type:"bar",
      x:sorted.map(d=>d.x), y:sorted.map(d=>d.y),
      error_y: errY(sorted),
      marker:{ color:colors[0], opacity:.85, line:{color:"black",width:0.8} },
      text:      ST.labels ? sorted.map(d => (+d.y.toFixed(2)).toString()) : [],
      textposition: ST.labels ? "outside" : "none",
    }];
    // Jitter point overlay (single invisible-box trace)
    if (ST.showPoints !== "none") {
      const pts = data.map(d => ({ x:String(d[xCol]), y:parseFloat(d[yCol]) }))
                      .filter(d => isFinite(d.y));
      if (pts.length) {
        traces.push({
          type:"box",
          x:pts.map(d=>d.x), y:pts.map(d=>d.y),
          boxpoints:"all", jitter:0.4, pointpos:0,
          whiskerwidth:0, boxmean:false,
          fillcolor:"rgba(0,0,0,0)", line:{color:"rgba(0,0,0,0)",width:0},
          marker:{ color:ST.blackPoints ? "#333333" : colors[0], size:5, opacity:0.7,
                   line:{color:"white",width:0.5} },
          showlegend:false, name:"__pts__",
        });
      }
    }
  }
  return {
    traces,
    layout:{
      barmode:ST.barMode,
      boxmode:"group",
      xaxis:{ title:xCol, tickangle:-25, automargin:true, gridcolor:gridColor(), autorange:ST.revX?"reversed":true },
      yaxis:{ title:yLabel(yCol,aggFn), gridcolor:gridColor(), zeroline:true, autorange:ST.revY?"reversed":true },
      legend:{ title:{text:colorCol||""} },
      bargap:.25, bargroupgap:.1,
    },
  };
}

// ── Scatter ───────────────────────────────────────────────────────────
function renderScatter(data = DATA) {
  const xCol     = document.getElementById("selX").value;
  const yCol     = document.getElementById("selY").value;
  const colorCol = document.getElementById("selColor").value || null;
  const sizeCol  = document.getElementById("selSize").value  || null;
  if (!xCol || !yCol) return;

  const colors = getColors();
  let traces;

  if (colorCol) {
    const groups = [...new Set(data.map(d => d[colorCol]))];
    traces = groups.map((g, gi) => {
      const gd = data.filter(d => d[colorCol] === g);
      const sizes = sizeCol ? gd.map(d => Math.sqrt(Math.abs(parseFloat(d[sizeCol])||1))*5+4) : 8;
      return {
        type:"scatter", mode:"markers", name:String(g),
        x:gd.map(d=>parseFloat(d[xCol])), y:gd.map(d=>parseFloat(d[yCol])),
        marker:{ color:colors[gi%colors.length], size:sizes, opacity:.78, line:{color:"white",width:ST.lineThick*0.3} },
        hovertemplate:`<b>${xCol}:</b> %{x}<br><b>${yCol}:</b> %{y}<br><b>${colorCol}:</b> ${g}<extra></extra>`,
      };
    });
  } else {
    const xs    = data.map(d => parseFloat(d[xCol]));
    const ys    = data.map(d => parseFloat(d[yCol]));
    const sizes = sizeCol ? data.map(d=>Math.sqrt(Math.abs(parseFloat(d[sizeCol])||1))*5+4) : 8;
    traces = [{
      type:"scatter", mode:"markers",
      x:xs, y:ys,
      marker:{
        color:sizeCol?data.map(d=>parseFloat(d[sizeCol])):colors[0],
        colorscale:sizeCol?getSeqScale():undefined,
        showscale:!!sizeCol,
        colorbar:sizeCol?{title:{text:sizeCol,font:{size:ST.fontSize-1}},thickness:12,len:.7}:undefined,
        size:sizes, opacity:.78, line:{color:"white",width:ST.lineThick*0.3},
      },
      hovertemplate:`<b>${xCol}:</b> %{x}<br><b>${yCol}:</b> %{y}<extra></extra>`,
    }];
  }
  return {
    traces,
    layout:{
      xaxis:{ title:xCol, gridcolor:gridColor(), autorange:ST.revX?"reversed":true },
      yaxis:{ title:yCol, gridcolor:gridColor(), autorange:ST.revY?"reversed":true },
      legend:{ title:{text:colorCol||""} },
    },
  };
}

// ── Line ──────────────────────────────────────────────────────────────
function renderLine(data = DATA) {
  const xCol     = document.getElementById("selX").value;
  const yCol     = document.getElementById("selY").value;
  const colorCol = document.getElementById("selColor").value || null;
  const aggFn    = document.getElementById("selAgg").value;
  if (!xCol || !yCol) return;

  const colors = getColors();
  let traces;

  if (colorCol) {
    const agg    = aggregate(data, xCol, yCol, colorCol, aggFn);
    const groups = [...new Set(agg.map(d=>d.color))];
    traces = groups.map((g,gi) => {
      const gd = agg.filter(d=>d.color===g).sort((a,b)=>String(a.x).localeCompare(String(b.x)));
      return {
        type:"scatter", mode:"lines+markers", name:String(g),
        x:gd.map(d=>d.x), y:gd.map(d=>d.y),
        line:{color:colors[gi%colors.length], width:ST.lineThick*0.7},
        marker:{color:colors[gi%colors.length], size:6},
      };
    });
  } else {
    const agg = aggregate(data, xCol, yCol, null, aggFn).sort((a,b)=>String(a.x).localeCompare(String(b.x)));
    traces = [{
      type:"scatter", mode:"lines+markers",
      x:agg.map(d=>d.x), y:agg.map(d=>d.y),
      line:{color:colors[0], width:ST.lineThick*0.8},
      marker:{color:colors[0], size:6},
    }];
  }
  return {
    traces,
    layout:{
      xaxis:{ title:xCol, automargin:true, gridcolor:gridColor(), autorange:ST.revX?"reversed":true },
      yaxis:{ title:yLabel(yCol,aggFn), gridcolor:gridColor(), autorange:ST.revY?"reversed":true },
      legend:{ title:{text:colorCol||""} },
    },
  };
}

// ── Histogram ─────────────────────────────────────────────────────────
function renderHistogram(data = DATA) {
  const xCol     = document.getElementById("selX").value;
  const colorCol = document.getElementById("selColor").value || null;
  const bins     = parseInt(document.getElementById("binSlider").value) || 20;
  if (!xCol) return;

  const colors = getColors();
  let traces;

  if (colorCol) {
    const groups = [...new Set(data.map(d=>d[colorCol]))];
    traces = groups.map((g,gi) => ({
      type:"histogram", name:String(g),
      x:data.filter(d=>d[colorCol]===g).map(d=>parseFloat(d[xCol])),
      nbinsx:bins,
      marker:{ color:colors[gi%colors.length], opacity:.72, line:{color:"white",width:ST.lineThick*0.3} },
    }));
  } else {
    traces = [{
      type:"histogram",
      x:data.map(d=>parseFloat(d[xCol])),
      nbinsx:bins,
      marker:{ color:colors[0], opacity:.82, line:{color:"white",width:ST.lineThick*0.4} },
    }];
  }
  return {
    traces,
    layout:{
      barmode:"overlay",
      xaxis:{ title:xCol, gridcolor:gridColor(), autorange:ST.revX?"reversed":true },
      yaxis:{ title:"Count", gridcolor:gridColor(), autorange:ST.revY?"reversed":true },
      legend:{ title:{text:colorCol||""} },
    },
  };
}

// ── Box ───────────────────────────────────────────────────────────────
function renderBox(data = DATA) {
  const xCol     = document.getElementById("selX").value;
  const yCol     = document.getElementById("selY").value;
  const colorCol = document.getElementById("selColor").value || null;
  if (!xCol || !yCol) return;

  let traces;
  let isDiagonalBox = false;

  const ptsMode = ST.showPoints === "none" ? false : ST.showPoints;
  const ptJit   = ST.showPoints === "none" ? 0 : 0.35;
  const ptPos   = 0;   // pointpos:0 → points sit ON the box, Prism-style

  if (colorCol) {
    const groups = [...new Set(data.map(d=>String(d[colorCol])))];

    // Detect diagonal: each color group appears at exactly one x-position
    const colorXSetsBox = new Map();
    data.forEach(d => {
      const c = String(d[colorCol]);
      if (!colorXSetsBox.has(c)) colorXSetsBox.set(c, new Set());
      colorXSetsBox.get(c).add(String(d[xCol]));
    });
    isDiagonalBox = groups.every(g => (colorXSetsBox.get(g)?.size || 0) <= 1);

    if (isDiagonalBox) {
      // Build one box trace per x-group, colored by its corresponding color group
      // Map x-value → color group (first match wins)
      const xToColor = new Map();
      data.forEach(d => {
        const xv = String(d[xCol]);
        if (!xToColor.has(xv)) xToColor.set(xv, String(d[colorCol]));
      });
      const xGroups = [...new Set(data.map(d => String(d[xCol])))];
      traces = xGroups.map((xv) => {
        const cg    = xToColor.get(xv) || xv;
        const gi    = groups.indexOf(cg);
        const gIdx  = gi >= 0 ? gi : 0;
        const color = getGroupColor(cg, gIdx);
        return {
          type:"box", name:xv,
          y: data.filter(d=>String(d[xCol])===xv).map(d=>parseFloat(d[yCol])),
          fillcolor: color,
          marker:{ color: pointColor(cg, gIdx), size:4, opacity:0.65 },
          line:{ color: "black", width:ST.lineThick*0.5 },
          boxmean:true,
          boxpoints: ptsMode,
          jitter:    ptJit,
          pointpos:  ptPos,
        };
      });
    } else {
      traces = groups.map((g,gi) => {
        const gd = data.filter(d=>String(d[colorCol])===g);
        const groupColor = getGroupColor(String(g), gi);
        return {
          type:"box", name:String(g),
          x:gd.map(d=>String(d[xCol])),
          y:gd.map(d=>parseFloat(d[yCol])),
          fillcolor: groupColor,
          marker:{ color: pointColor(String(g), gi), size:4, opacity:0.65 },
          line:{ color: "black", width:ST.lineThick*0.5 },
          boxmean:true,
          boxpoints: ptsMode,
          jitter:    ptJit,
          pointpos:  ptPos,
        };
      });
    }
  } else {
    const groups = [...new Set(data.map(d=>d[xCol]))];
    traces = groups.map((g,gi) => {
      const groupColor = getGroupColor(String(g), gi);
      return {
        type:"box", name:String(g),
        y:data.filter(d=>d[xCol]===g).map(d=>parseFloat(d[yCol])),
        fillcolor: groupColor,
        marker:{ color: pointColor(String(g), gi), size:4, opacity:0.65 },
        line:{ color: "black", width:ST.lineThick*0.5 },
        boxmean:true,
        boxpoints: ptsMode,
        jitter:    ptJit,
        pointpos:  ptPos,
      };
    });
  }
  return {
    traces,
    layout:{
      // Diagonal (xCol===colorCol): "overlay" centers each box on its x-slot
      // Grouped: "group" places colored boxes side by side per x-category
      boxmode: isDiagonalBox ? "overlay" : "group",
      xaxis:{ title:xCol, automargin:true, gridcolor:gridColor(), autorange:ST.revX?"reversed":true },
      yaxis:{ title:yCol, gridcolor:gridColor(), autorange:ST.revY?"reversed":true },
      legend:{ title:{text:colorCol||""} },
    },
  };
}

// ── Violin ───────────────────────────────────────────────────────────────
function renderViolin(data = DATA) {
  const xCol     = document.getElementById("selX").value;
  const yCol     = document.getElementById("selY").value;
  const colorCol = document.getElementById("selColor").value || null;
  if (!xCol || !yCol) return;

  const colors  = getColors();
  const ptsMode = ST.showPoints === "none" ? false : ST.showPoints;
  const ptJit   = ST.showPoints === "none" ? 0 : 0.35;
  const ptPos   = 0;   // center points on the violin, Prism-style

  let traces;
  let isDiagonalViolin = false;

  if (colorCol) {
    const groups = [...new Set(data.map(d=>d[colorCol]))];

    // Detect diagonal: each color group appears at exactly one x-position
    const colorXSetsVio = new Map();
    data.forEach(d => {
      const c = String(d[colorCol]);
      if (!colorXSetsVio.has(c)) colorXSetsVio.set(c, new Set());
      colorXSetsVio.get(c).add(String(d[xCol]));
    });
    isDiagonalViolin = groups.every(g => (colorXSetsVio.get(g)?.size || 0) <= 1);

    if (isDiagonalViolin) {
      // One violin per x-group, y-only, colored by its color group → violinmode:"overlay" centers it
      const xToColor = new Map();
      data.forEach(d => {
        const xv = String(d[xCol]);
        if (!xToColor.has(xv)) xToColor.set(xv, String(d[colorCol]));
      });
      const xGroups = [...new Set(data.map(d => String(d[xCol])))];
      traces = xGroups.map(xv => {
        const cg   = xToColor.get(xv) || xv;
        const gi   = groups.indexOf(cg);
        const gIdx = gi >= 0 ? gi : 0;
        const c    = getGroupColor(cg, gIdx);
        return {
          type: "violin", name: xv,
          y: data.filter(d=>String(d[xCol])===xv).map(d=>parseFloat(d[yCol])),
          fillcolor: hexToRgba(c, 0.45),
          line:    { color: "black", width: ST.lineThick * 0.5 },
          box:     { visible: true },
          meanline:{ visible: true },
          points:  ptsMode,
          jitter:  ptJit,
          pointpos: ptPos,
          marker:  { color: pointColor(cg, gIdx), size: 3, opacity: 0.65 },
        };
      });
    } else {
      traces = groups.map((g, gi) => {
        const gd = data.filter(d=>d[colorCol]===g);
        const c  = getGroupColor(String(g), gi);
        return {
          type: "violin", name: String(g),
          x: gd.map(d=>String(d[xCol])),
          y: gd.map(d=>parseFloat(d[yCol])),
          fillcolor: hexToRgba(c, 0.45),
          line:    { color: "black", width: ST.lineThick * 0.5 },
          box:     { visible: true },
          meanline:{ visible: true },
          points:  ptsMode,
          jitter:  ptJit,
          pointpos: ptPos,
          marker:  { color: pointColor(String(g), gi), size: 3, opacity: 0.65 },
        };
      });
    }
  } else {
    const groups = [...new Set(data.map(d=>d[xCol]))];
    traces = groups.map((g, gi) => {
      const c = getGroupColor(String(g), gi);
      return {
        type: "violin", name: String(g),
        y: data.filter(d=>d[xCol]===g).map(d=>parseFloat(d[yCol])),
        fillcolor: hexToRgba(c, 0.45),
        line:    { color: "black", width: ST.lineThick * 0.5 },
        box:     { visible: true },
        meanline:{ visible: true },
        points:  ptsMode,
        jitter:  ptJit,
        pointpos: ptPos,
        marker:  { color: pointColor(String(g), gi), size: 3, opacity: 0.65 },
      };
    });
  }
  return {
    traces,
    layout: {
      // Diagonal (xCol===colorCol): "overlay" centers each violin on its x-slot
      violinmode: isDiagonalViolin ? "overlay" : "group",
      xaxis: { title: xCol, automargin: true, gridcolor: gridColor(), autorange: ST.revX ? "reversed" : true },
      yaxis: { title: yCol, gridcolor: gridColor(), autorange: ST.revY ? "reversed" : true },
      legend:{ title: { text: colorCol || "" } },
    },
  };
}

// ── Heatmap — long format (single or faceted) ─────────────────────────
function renderHeatmap(data = DATA) {
  const xCol     = document.getElementById("selX").value;
  const yCol     = document.getElementById("selY").value;
  const zCol     = document.getElementById("selColor").value;
  const facetCol = document.getElementById("selFacet").value || null;
  const aggFn    = document.getElementById("selAgg").value;
  if (!xCol || !yCol || !zCol) return;

  const cbConf = (xPos) => ({
    title:{ text:zCol, font:{size:ST.fontSize-1} },
    x:xPos, thickness:13, tickfont:{size:ST.fontSize-2}, len:0.85,
  });

  if (!facetCol) {
    const { xVals, yVals, z } = buildHeatMatrix(data, xCol, yCol, zCol, aggFn);
    return {
      traces:[{
        type:"heatmap", x:xVals, y:yVals, z,
        colorscale:getHeatScale(), colorbar:cbConf(1.02),
        hoverongaps:false,
        hovertemplate:`${xCol}: %{x}<br>${yCol}: %{y}<br>${zCol}: %{z:.3g}<extra></extra>`,
      }],
      layout:{
        margin:{ l:110, r:90, t:30, b:80 },
        xaxis:{ title:xCol, tickangle:-30, automargin:true, showgrid:false, autorange:ST.revX?"reversed":true },
        yaxis:{ title:yCol, automargin:true, showgrid:false, autorange:ST.revY?"reversed":true },
      },
    };
  }

  // Faceted heatmap: use renderFaceted wrapper
  // (We strip facetCol from selFacet temporarily to avoid infinite recursion)
  const nCols = parseInt(document.getElementById("facetNCols")?.value) || ST.nFacetCols;
  const allZ  = data.map(d=>parseFloat(d[zCol])).filter(v=>!isNaN(v));
  const zMin  = Math.min(...allZ), zMax = Math.max(...allZ);

  return renderFaceted(subset => {
    const { xVals, yVals, z } = buildHeatMatrix(subset, xCol, yCol, zCol, aggFn);
    return {
      traces:[{
        type:"heatmap", x:xVals, y:yVals, z,
        colorscale:getHeatScale(), zmin:zMin, zmax:zMax,
        showscale:true, colorbar:cbConf(1.01),
        hoverongaps:false,
        hovertemplate:`${xCol}: %{x}<br>${yCol}: %{y}<br>${zCol}: %{z:.3g}<extra></extra>`,
      }],
      layout:{
        margin:{ l:110, r:90, t:30, b:80 },
        xaxis:{ tickangle:-30, automargin:true, showgrid:false },
        yaxis:{ automargin:true, showgrid:false },
      },
    };
  }, data, facetCol, nCols, false);
}

// ── Wide format heatmap helpers ───────────────────────────────────────
function makeDiscreteColorscale(n, colors) {
  if (n <= 1) return [[0, colors[0]||"#ccc"],[1, colors[0]||"#ccc"]];
  const steps = [];
  for (let i = 0; i < n; i++) {
    steps.push([i / n, colors[i % colors.length]]);
    if (i < n - 1) steps.push([(i + 1) / n - 1e-9, colors[i % colors.length]]);
  }
  steps.push([1, colors[(n - 1) % colors.length]]);
  return steps;
}

// ── Heatmap — wide format ─────────────────────────────────────────────
function renderWideHeatmap() {
  // Row label: use first meta col, or fall back to row index string
  const rowLabelCol = document.getElementById("selRowLabel")?.value || WS.metaCols[0] || null;
  const rowSplitCol = document.getElementById("selRowSplit")?.value || null;
  const doLog2      = ST.log2;
  const doCenter    = ST.center;
  const doColSplit  = ST.colSplit && WS.groupOrder.length > 1;

  if (!WS.sampleCols.length) return;

  // 1. Sort rows by rowSplitCol (stable sort)
  let rows = DATA.map((r, i) => ({ r, i }));
  if (rowSplitCol) {
    rows.sort((a, b) =>
      String(a.r[rowSplitCol]).localeCompare(String(b.r[rowSplitCol])) || a.i - b.i
    );
  }
  rows = rows.map(x => x.r);

  // Apply row filter (gene/feature filter); null = no filter, empty Set = show nothing
  if (WS.rowFilter !== null) {
    rows = rows.filter(r => {
      const label = rowLabelCol ? String(r[rowLabelCol] ?? "") : "";
      return WS.rowFilter.has(label);
    });
  }

  // 2. Feature names — use label col or row index
  const featureNames = rows.map((r, i) =>
    rowLabelCol ? String(r[rowLabelCol] ?? "") : String(i + 1)
  );

  // 3. Sample columns in group order
  let sampleCols;
  if (doColSplit) {
    sampleCols = WS.groupOrder.flatMap(g =>
      WS.sampleCols.filter(s => WS.sampleGroups[s] === g)
    );
  } else {
    sampleCols = [...WS.sampleCols];
  }
  // Apply sample filter — null = no filter, empty Set = show nothing
  if (WS.sampleFilter !== null) {
    sampleCols = sampleCols.filter(sc => WS.sampleFilter.has(sc));
  }

  const N = featureNames.length;
  const M = sampleCols.length;
  if (!N || !M) return;

  // 4. Build z matrix
  let z = rows.map(row =>
    sampleCols.map(sc => {
      const v = parseFloat(row[sc]);
      return isNaN(v) ? null : v;
    })
  );

  // 5. Transforms
  if (doLog2)  z = z.map(row => row.map(v => v === null ? null : Math.log2(Math.max(0, v) + 1)));
  if (doCenter) {
    z = z.map(row => {
      const vals = row.filter(v => v !== null);
      if (!vals.length) return row;
      const mean = vals.reduce((a,b) => a + b, 0) / vals.length;
      return row.map(v => v === null ? null : v - mean);
    });
  }

  // 6. Color scale
  const allZ = z.flat().filter(v => v !== null);
  if (!allZ.length) return;
  let zmin, zmax;
  if (doCenter) {
    // Symmetric range so the scale midpoint sits at 0
    const maxAbs = Math.max(...allZ.map(Math.abs));
    zmin = -maxAbs; zmax = maxAbs;
  } else {
    zmin = Math.min(...allZ); zmax = Math.max(...allZ);
  }
  // Colorscale follows the user's palette choice (Center auto-switches it to "div")
  const colorscale = getHeatScale();

  // 7. Layout fractions
  const hasRowAnn = !!rowSplitCol;
  const hasColAnn = doColSplit;
  const TANN_H    = hasColAnn ? 0.04 : 0;
  const GAP       = 0.012;
  // Shrink main plot when there are row annotations to leave room for labels + colorbar
  const MAIN_X1   = hasRowAnn ? 0.70 : 0.88;
  const RANN_W    = hasRowAnn ? 0.04 : 0;
  const RANN_X0   = hasRowAnn ? MAIN_X1 + GAP : 0;
  const RANN_X1   = hasRowAnn ? RANN_X0 + RANN_W : 0;
  const MAIN_Y0   = 0;
  const MAIN_Y1   = hasColAnn ? 1 - TANN_H - GAP : 1.0;

  // Colorbar label
  const userLabel = (document.getElementById("valLabel")?.value || "").trim();
  const autoLabel = doCenter ? "Centered" : (doLog2 ? "log₂(x+1)" : "Value");
  const cbLabel   = userLabel || autoLabel;

  // Compute max row-group label length for spacing
  let maxLabelChars = 0;
  if (hasRowAnn) {
    const rowGroups = [...new Set(rows.map(r => String(r[rowSplitCol])))];
    maxLabelChars   = Math.max(...rowGroups.map(g => g.length), 1);
  }

  // Colorbar x when right:
  // Labels start at RANN_X1 + 0.01 (xanchor left) and span ~0.0085 paper units per char.
  // Place colorbar past the end of labels with an extra gap.
  const LABEL_X0      = hasRowAnn ? RANN_X1 + 0.01 : 0;
  const CHAR_W_PAPER  = 0.0085;   // paper fractions per character
  const CB_AFTER_GAP  = 0.020;    // gap between label end and colorbar
  const cbXRight      = hasRowAnn
    ? LABEL_X0 + maxLabelChars * CHAR_W_PAPER + CB_AFTER_GAP
    : MAIN_X1 + 0.025;

  // Colorbar configuration
  const cbRight = ST.cbPosition !== "bottom";
  const cbConf  = cbRight
    ? {
        title:     { text: cbLabel, font:{ size: Math.max(8, ST.fontSize - 2) } },
        x:         cbXRight,
        xanchor:   "left",
        thickness: 14,
        tickfont:  { size: Math.max(7, ST.fontSize - 3) },
        len:       0.70,
      }
    : {
        orientation: "h",
        title:       { text: cbLabel, font:{ size: Math.max(8, ST.fontSize - 2) }, side:"bottom" },
        x:           0.5, xanchor: "center",
        y:           -0.18, yanchor: "top",
        thickness:   14,
        tickfont:    { size: Math.max(7, ST.fontSize - 3) },
        len:         0.60,
      };

  const traces      = [];
  const annotations = [];
  const shapes      = [];
  const layoutAxes  = {};
  const CAT_COLORS  = getColors();

  // 8. Main heatmap
  const tickFontSz = Math.max(7, Math.min(ST.fontSize, Math.floor(600 / N)));
  traces.push({
    type: "heatmap",
    x: sampleCols, y: featureNames, z,
    colorscale, zmin, zmax,
    xaxis: "x", yaxis: "y",
    colorbar: cbConf,
    hoverongaps: false,
    hovertemplate: `%{y}<br>%{x}: %{z:.3g}<extra></extra>`,
  });

  layoutAxes.xaxis = {
    domain: [0, MAIN_X1],
    showticklabels: M <= 30,
    tickangle: -45, automargin: true, showgrid: false,
    tickfont: { size: ST.fontSize - 2 },
  };
  layoutAxes.yaxis = {
    domain: [MAIN_Y0, MAIN_Y1],
    automargin: true, showgrid: false,
    tickfont: { size: tickFontSz },
    autorange: "reversed",
  };

  function buildSegments(arr, keyFn) {
    const segs = [];
    arr.forEach((item, i) => {
      const k = keyFn(item);
      if (!segs.length || segs[segs.length - 1].key !== k)
        segs.push({ key: k, start: i, end: i + 1 });
      else
        segs[segs.length - 1].end = i + 1;
    });
    return segs;
  }

  // 9. Row annotation strip
  if (hasRowAnn) {
    const rowSegs   = buildSegments(rows, r => String(r[rowSplitCol]));
    const rowGroups = [...new Set(rowSegs.map(s => s.key))];
    const nRG       = rowGroups.length;
    const rowIdx    = rows.map(r => rowGroups.indexOf(String(r[rowSplitCol])));

    traces.push({
      type: "heatmap",
      x: [rowSplitCol], y: featureNames,
      z: rowIdx.map(i => [i]),
      text: rows.map(r => [String(r[rowSplitCol])]),
      zmin: -0.5, zmax: nRG - 0.5,
      colorscale: makeDiscreteColorscale(nRG, CAT_COLORS),
      showscale: false,
      xaxis: "x2", yaxis: "y",
      hovertemplate: `%{y}<br>${rowSplitCol}: %{text}<extra></extra>`,
    });

    layoutAxes.xaxis2 = {
      domain: [RANN_X0, RANN_X1],
      showticklabels: false, showgrid: false, anchor: "y",
    };

    rowSegs.forEach(seg => {
      const gi   = rowGroups.indexOf(seg.key);
      const topY = MAIN_Y1 - (MAIN_Y1 - MAIN_Y0) * seg.start / N;
      const botY = MAIN_Y1 - (MAIN_Y1 - MAIN_Y0) * seg.end   / N;
      const midY = (topY + botY) / 2;

      annotations.push({
        text: `<b>${seg.key}</b>`,
        xref:"paper", yref:"paper",
        x: RANN_X1 + 0.01, y: midY,
        xanchor:"left", yanchor:"middle",
        showarrow: false,
        font: { size: Math.max(7, ST.fontSize - 3), color:CAT_COLORS[gi % CAT_COLORS.length], family:FONT },
      });

      if (seg.end < N) {
        shapes.push({
          type:"line", xref:"paper", yref:"paper",
          x0:0, x1: RANN_X1 + 0.002,
          y0:botY, y1:botY,
          line:{ color:"white", width: ST.lineThick }, layer:"above",
        });
      }
    });
  }

  // 10. Top annotation strip
  if (hasColAnn) {
    const colSegs    = buildSegments(sampleCols, s => WS.sampleGroups[s]);
    const nCG        = WS.groupOrder.length;
    const colIdx     = sampleCols.map(s => WS.groupOrder.indexOf(WS.sampleGroups[s]));

    traces.push({
      type: "heatmap",
      x: sampleCols, y: [""],
      z: [colIdx],
      text: [sampleCols.map(s => WS.sampleGroups[s])],
      zmin: -0.5, zmax: nCG - 0.5,
      colorscale: makeDiscreteColorscale(nCG, CAT_COLORS),
      showscale: false,
      xaxis: "x", yaxis: "y2",
      hovertemplate: `%{x}<br>Group: %{text}<extra></extra>`,
    });

    layoutAxes.yaxis2 = {
      domain: [MAIN_Y1 + GAP, 1.0],
      showticklabels: false, showgrid: false, anchor: "x",
    };

    colSegs.forEach(seg => {
      const gi    = WS.groupOrder.indexOf(seg.key);
      const leftX = (MAIN_X1) * seg.start / M;
      const ritX  = (MAIN_X1) * seg.end   / M;
      const midX  = (leftX + ritX) / 2;

      annotations.push({
        text: `<b>${seg.key}</b>`,
        xref:"paper", yref:"paper",
        x: midX, y: 1.01,
        xanchor:"center", yanchor:"bottom",
        showarrow: false,
        font: { size: Math.max(8, ST.fontSize - 2), color:CAT_COLORS[gi % CAT_COLORS.length], family:FONT },
      });

      if (seg.end < M) {
        shapes.push({
          type:"line", xref:"paper", yref:"paper",
          x0:ritX, x1:ritX,
          y0: MAIN_Y0 - 0.01, y1: 1.02,
          line:{ color:"white", width: ST.lineThick }, layer:"above",
        });
      }
    });
  }

  // Right margin: must fit row labels + colorbar (or just colorbar if no row ann)
  // Labels span maxLabelChars * ~8px + 20px padding
  const labelPx   = hasRowAnn ? (maxLabelChars * 8 + 24) : 0;
  const cbPx      = cbRight ? 80 : 0;   // colorbar + ticks + title
  const rMargin   = Math.max(60, labelPx + cbPx + 20);

  // Bottom margin: needs room for bottom colorbar
  const bMargin   = cbRight ? (M <= 30 ? 80 : 20) : 150;

  // ── Column annotation tracks (sample metadata fields) ──────────────
  const annotFields = (WS.metaFields || []).slice(0, 5); // max 5 tracks
  if (annotFields.length && Object.keys(WS.sampleMeta || {}).length) {
    const TRACK_H = 0.045;
    const GAP_T   = 0.010;
    const totalAnnotH = annotFields.length * (TRACK_H + GAP_T);
    const mainTop = Math.max(0.4, 1 - totalAnnotH - 0.02);

    // Shrink main heatmap yaxis domain to leave room at top for annotation tracks
    layoutAxes.yaxis = { ...(layoutAxes.yaxis || {}), domain: [MAIN_Y0, mainTop] };
    // Shift existing column annotation strip (yaxis2) if present
    if (hasColAnn && layoutAxes.yaxis2) {
      layoutAxes.yaxis2 = { ...layoutAxes.yaxis2,
        domain: [mainTop + GAP, mainTop + GAP + TANN_H] };
    }

    const annotColors = getColors();
    annotFields.forEach((field, fi) => {
      // Each annotation track sits above the main heatmap
      const trackBottom = 1 - (fi + 1) * (TRACK_H + GAP_T) + GAP_T;
      const trackTop    = trackBottom + TRACK_H;
      const axIdx = fi + (hasColAnn ? 3 : 2); // avoid collision with yaxis2 (col split)
      const yKey  = `yaxis${axIdx}`;
      const yRef  = `y${axIdx}`;

      const uniqueVals = [...new Set(sampleCols.map(sc => String(WS.sampleMeta[sc]?.[field] ?? "")))];
      const n = uniqueVals.length;
      // Build a stepped colorscale
      const cs = [];
      uniqueVals.forEach((_, i) => {
        const c = annotColors[i % annotColors.length];
        if (i === 0) cs.push([0, c]);
        cs.push([Math.min(1, Math.max(0, (i + 0.5) / Math.max(1, n))), c]);
        if (i === n - 1) cs.push([1, c]);
      });
      const colorscaleAnn = cs.length >= 2 ? cs : [[0, annotColors[0]], [1, annotColors[0]]];
      const zRow = sampleCols.map(sc => uniqueVals.indexOf(String(WS.sampleMeta[sc]?.[field] ?? "")));

      traces.push({
        type: "heatmap",
        x: sampleCols,
        y: [field],
        z: [zRow],
        colorscale: colorscaleAnn,
        showscale: false,
        xaxis: "x",
        yaxis: yRef,
        hovertemplate: `%{x}<br>${esc(field)}: %{text}<extra></extra>`,
        text: [sampleCols.map(sc => String(WS.sampleMeta[sc]?.[field] ?? "—"))],
        zmin: 0, zmax: Math.max(1, n - 1),
      });

      layoutAxes[yKey] = {
        domain: [trackBottom, trackTop],
        anchor: "x",
        showticklabels: true,
        tickfont: { size: 9, color: "#555" },
        showgrid: false,
        zeroline: false,
        fixedrange: true,
      };
    });
  }

  return {
    traces,
    layout: {
      annotations, shapes,
      margin: {
        l: 140,
        r: rMargin,
        t: annotFields.length ? Math.max(hasColAnn ? 60 : 30, 20 + annotFields.length * 5) : (hasColAnn ? 60 : 30),
        b: bMargin,
      },
      paper_bgcolor:"white", plot_bgcolor:"white",
      ...layoutAxes,
    },
  };
}

// ═══════════════════════════════════════════════════════════════════
//  ALLUVIAL / SANKEY
//  X column = Source (left layer)   Y column = Target (right layer)
//  Color column = optional Middle node → creates 3-layer sankey
//  Size column  = optional value/weight (defaults to count = 1 per row)
// ═══════════════════════════════════════════════════════════════════
function renderAlluvial(data = DATA) {
  const xCol   = document.getElementById("selX").value;
  const yCol   = document.getElementById("selY").value;
  const midCol = document.getElementById("selColor").value || null;
  const valCol = document.getElementById("selSize").value  || null;

  if (!xCol || !yCol || xCol === yCol) return null;

  const colors = getColors();
  const has3   = !!(midCol && midCol !== xCol && midCol !== yCol);

  // Unique node values per layer (stable sort)
  const srcVals = [...new Set(data.map(r => String(r[xCol] ?? "")))].filter(Boolean).sort();
  const tgtVals = [...new Set(data.map(r => String(r[yCol] ?? "")))].filter(Boolean).sort();
  const midVals = has3 ? [...new Set(data.map(r => String(r[midCol] ?? "")))].filter(Boolean).sort() : [];

  if (!srcVals.length || !tgtVals.length) return null;

  const nS = srcVals.length, nM = midVals.length;
  const srcIdx = v => srcVals.indexOf(v);
  const midIdx = v => nS + midVals.indexOf(v);
  const tgtIdx = v => nS + nM + tgtVals.indexOf(v);

  const srcColors = srcVals.map((_, i) => colors[i % colors.length]);
  const midColors = midVals.map((_, i) => colors[i % colors.length]);
  const tgtColors = tgtVals.map((_, i) => colors[(nS + i) % colors.length]);

  const nodeLabels = [...srcVals, ...midVals, ...tgtVals];
  const nodeColors = [...srcColors, ...midColors, ...tgtColors];

  const lnkSrc = [], lnkTgt = [], lnkVal = [], lnkCol = [], lnkLbl = [];

  if (has3) {
    // 3-layer: src → mid, then mid → tgt (aggregate each leg separately)
    const L1 = new Map(), L2 = new Map();
    data.forEach(r => {
      const sv = String(r[xCol] ?? ""), mv = String(r[midCol] ?? ""), tv = String(r[yCol] ?? "");
      if (!sv || !mv || !tv) return;
      const v = valCol ? (parseFloat(r[valCol]) || 0) : 1;
      L1.set(`${sv}\x00${mv}`, (L1.get(`${sv}\x00${mv}`) || 0) + v);
      L2.set(`${mv}\x00${tv}`, (L2.get(`${mv}\x00${tv}`) || 0) + v);
    });
    L1.forEach((v, key) => {
      const [sv, mv] = key.split("\x00");
      const si = srcIdx(sv), mi = midIdx(mv);
      if (si < 0 || mi < 0) return;
      lnkSrc.push(si); lnkTgt.push(mi); lnkVal.push(v);
      lnkCol.push(hexToRgba(srcColors[si], 0.40));
      lnkLbl.push(`${sv} → ${mv}`);
    });
    L2.forEach((v, key) => {
      const [mv, tv] = key.split("\x00");
      const mi = midIdx(mv), ti = tgtIdx(tv);
      if (mi < 0 || ti < 0) return;
      lnkSrc.push(mi); lnkTgt.push(ti); lnkVal.push(v);
      lnkCol.push(hexToRgba(midColors[midVals.indexOf(mv)], 0.40));
      lnkLbl.push(`${mv} → ${tv}`);
    });
  } else {
    // 2-layer: src → tgt
    const L = new Map();
    data.forEach(r => {
      const sv = String(r[xCol] ?? ""), tv = String(r[yCol] ?? "");
      if (!sv || !tv) return;
      const v = valCol ? (parseFloat(r[valCol]) || 0) : 1;
      L.set(`${sv}\x00${tv}`, (L.get(`${sv}\x00${tv}`) || 0) + v);
    });
    L.forEach((v, key) => {
      const [sv, tv] = key.split("\x00");
      const si = srcIdx(sv), ti = tgtIdx(tv);
      if (si < 0 || ti < 0) return;
      lnkSrc.push(si); lnkTgt.push(ti); lnkVal.push(v);
      lnkCol.push(hexToRgba(srcColors[si], 0.40));
      lnkLbl.push(`${sv} → ${tv}`);
    });
  }

  if (!lnkVal.length) return null;

  const subtitle = has3
    ? `← <b>${xCol}</b> &nbsp;&nbsp; ${midCol} &nbsp;&nbsp; <b>${yCol}</b> →`
    : `← <b>${xCol}</b> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <b>${yCol}</b> →`;

  return {
    traces: [{
      type: "sankey", orientation: "h",
      node: {
        pad: 15, thickness: 20,
        line: { color:"white", width:0.5 },
        label: nodeLabels, color: nodeColors,
        hovertemplate: "%{label}<br>Total flow: %{value}<extra></extra>",
      },
      link: {
        source: lnkSrc, target: lnkTgt, value: lnkVal,
        color:  lnkCol, label: lnkLbl,
        hovertemplate: "%{label}<br>Value: %{value:.3g}<extra></extra>",
      },
    }],
    layout: {
      margin: { l:30, r:30, t:52, b:20 },
      annotations: [{
        text: subtitle,
        xref:"paper", yref:"paper", x:0.5, y:1.06,
        showarrow:false, xanchor:"center",
        font:{ size: ST.fontSize, color:"#7a8594" },
      }],
    },
  };
}

// ═══════════════════════════════════════════════════════════════════
//  EXPORT  (Problem 7: uses ST.exportW / ST.exportH)
// ═══════════════════════════════════════════════════════════════════
async function exportFig(fmt) {
  const W = ST.exportW || 2400;
  const H = ST.exportH || 1600;

  if (fmt === "pdf") {
    showSpin("Generating PDF…");
    try {
      const svgUrl = await Plotly.toImage("plot", { format:"svg", width:W, height:H });
      const svgStr = decodeURIComponent(svgUrl.split(",").slice(1).join(","));
      const svgEl  = new DOMParser().parseFromString(svgStr, "image/svg+xml").documentElement;
      svgEl.setAttribute("width", W); svgEl.setAttribute("height", H);
      const { jsPDF } = window.jspdf;
      const pdf = new jsPDF({ orientation: W>H ? "landscape":"portrait", unit:"pt", format:[W, H] });
      await svg2pdf(svgEl, pdf, { x:0, y:0, width:W, height:H });
      pdf.save("chart.pdf");
    } catch(e) {
      const png = await Plotly.toImage("plot", { format:"png", width:W, height:H, scale:2 });
      const { jsPDF } = window.jspdf;
      const pdf = new jsPDF({ orientation: W>H ? "landscape":"portrait", unit:"pt", format:"a3" });
      const pw = pdf.internal.pageSize.getWidth();
      const ph = pdf.internal.pageSize.getHeight();
      pdf.addImage(png, "PNG", 0, 0, pw, ph);
      pdf.save("chart.pdf");
    }
    hideSpin(); return;
  }
  // scale:1 keeps fonts proportional to the chart — scale:2 would double canvas size
  // making fonts appear half as large relative to the figure
  Plotly.downloadImage("plot", {
    format: fmt, filename:"chart",
    width:  W,
    height: H,
    scale:  1,
  });
}
