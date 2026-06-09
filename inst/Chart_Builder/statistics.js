/*! © 2024 Nourhan Abdelfattah — scSidekick R package
 *  Bundled for use via scSidekick::ChartBuilder() only.
 *  Not licensed for standalone redistribution.
 */
// ═══════════════════════════════════════════════════════════════════
//  DATA HELPERS  (depends on: config.js)
// ═══════════════════════════════════════════════════════════════════

function aggregate(rows, xCol, yCol, colorCol, aggFn) {
  const groups = new Map();
  rows.forEach(row => {
    const xv  = String(row[xCol] ?? "(empty)");
    const cv  = colorCol ? String(row[colorCol] ?? "(empty)") : "__all__";
    const key = `${xv}\x00${cv}`;
    if (!groups.has(key)) groups.set(key, { x:xv, color:cv, vals:[], count:0 });
    const g = groups.get(key);
    g.count++;
    if (yCol) { const v = parseFloat(row[yCol]); if (!isNaN(v)) g.vals.push(v); }
  });
  return [...groups.values()].map(g => {
    let y, sd = 0, sem = 0;
    const v = g.vals, n = v.length;
    switch (aggFn) {
      case "sum":   y = v.reduce((a,b)=>a+b,0); break;
      case "mean":
        y = n ? v.reduce((a,b)=>a+b,0)/n : 0;
        if (n > 1) {
          sd  = Math.sqrt(v.reduce((a,b)=>a+(b-y)**2,0)/(n-1));
          sem = sd / Math.sqrt(n);
        }
        break;
      case "count": y = g.count; break;
      case "max":   y = n ? Math.max(...v) : 0; break;
      case "min":   y = n ? Math.min(...v) : 0; break;
      default:      y = v[0] ?? g.count;
    }
    return {
      x: g.x, y,
      color: g.color === "__all__" ? null : g.color,
      n, sd: +sd.toFixed(6), sem: +sem.toFixed(6),
    };
  });
}

function sortAgg(agg, mode) {
  if (mode === "alpha") return [...agg].sort((a,b) => String(a.x).localeCompare(String(b.x)));
  if (mode === "asc")   return [...agg].sort((a,b) => a.y - b.y);
  if (mode === "desc")  return [...agg].sort((a,b) => b.y - a.y);
  return agg;
}

function uniqueXOrder(agg, mode) {
  const xVals = [...new Set(agg.map(d=>d.x))];
  if (mode === "alpha") return xVals.sort((a,b)=>String(a).localeCompare(String(b)));
  if (mode === "asc" || mode === "desc") {
    const tot = {};
    agg.forEach(d => { tot[d.x] = (tot[d.x]||0) + d.y; });
    return xVals.sort((a,b) => mode==="asc" ? (tot[a]||0)-(tot[b]||0) : (tot[b]||0)-(tot[a]||0));
  }
  return xVals;
}

function buildHeatMatrix(rows, xCol, yCol, zCol, aggFn) {
  const xVals = [...new Set(rows.map(r => String(r[xCol]||"(empty)")))];
  const yVals = [...new Set(rows.map(r => String(r[yCol]||"(empty)")))];
  const lk = {};
  rows.forEach(r => {
    const k = `${r[xCol]}\x00${r[yCol]}`;
    if (!lk[k]) lk[k] = [];
    const v = parseFloat(r[zCol]); if (!isNaN(v)) lk[k].push(v);
  });
  const z = yVals.map(yv => xVals.map(xv => {
    const vs = lk[`${xv}\x00${yv}`] || [];
    if (!vs.length) return null;
    switch(aggFn) {
      case "sum":  return vs.reduce((a,b)=>a+b,0);
      case "mean": return vs.reduce((a,b)=>a+b,0)/vs.length;
      case "max":  return Math.max(...vs);
      case "min":  return Math.min(...vs);
      default:     return vs[0];
    }
  }));
  return { xVals, yVals, z };
}

// ═══════════════════════════════════════════════════════════════════
//  WIDE → LONG CONVERSION
// ═══════════════════════════════════════════════════════════════════

// Pivot wide-format DATA into a long (tidy) format for non-heatmap renderers.
// Produces columns: [all WS.metaCols..., "Sample", "Group" (if >1 group), "Value"]
function wideToLong() {
  const hasGroups  = WS.groupOrder.length > 1;
  const rowLabelCol = WS.metaCols[0] || null;   // primary row-label column for rowFilter
  const rows = [];
  DATA.forEach(row => {
    // Apply row filter — null = no filter, empty Set = show nothing
    if (WS.rowFilter !== null && rowLabelCol) {
      const label = String(row[rowLabelCol] ?? "");
      if (!WS.rowFilter.has(label)) return;
    }
    const meta = {};
    WS.metaCols.forEach(mc => { meta[mc] = row[mc]; });
    WS.sampleCols.forEach(sc => {
      // Apply sampleFilter — null = no filter, empty Set = show nothing
      if (WS.sampleFilter !== null && !WS.sampleFilter.has(sc)) return;
      const v = parseFloat(row[sc]);
      const longRow = {
        ...meta,
        Sample: sc,
        Value:  isNaN(v) ? null : v,
      };
      if (hasGroups) longRow.Group = WS.sampleGroups[sc] || sc;
      // Include sample metadata fields
      if (WS.metaFields?.length && WS.sampleMeta?.[sc]) {
        WS.metaFields.forEach(f => { longRow[f] = WS.sampleMeta[sc][f]; });
      }
      rows.push(longRow);
    });
  });
  return rows;
}

// Returns a COLS-compatible array describing the melted data structure.
function getMeltedCols() {
  const cols = [];
  WS.metaCols.forEach(name => {
    const ci = COLS.find(c => c.name === name);
    cols.push(ci || { name, type:"cat", uniq:1, min:null, max:null });
  });
  cols.push({ name:"Sample", type:"cat", uniq:WS.sampleCols.length, min:null, max:null });
  if (WS.groupOrder.length > 1) {
    cols.push({ name:"Group", type:"cat", uniq:WS.groupOrder.length, min:null, max:null });
  }
  const vals = DATA.flatMap(r => WS.sampleCols.map(sc => parseFloat(r[sc]))).filter(v => !isNaN(v));
  const vMin = vals.length ? Math.min(...vals) : null;
  const vMax = vals.length ? Math.max(...vals) : null;
  cols.push({ name:"Value", type:"num", uniq:new Set(vals).size, min:vMin, max:vMax });
  // Include sample metadata fields
  if (WS.metaFields?.length) {
    WS.metaFields.forEach(f => {
      const fVals = WS.sampleCols.map(sc => WS.sampleMeta?.[sc]?.[f]).filter(v => v != null);
      const uniq = [...new Set(fVals)];
      cols.push({ name:f, type:"cat", uniq:uniq.length, min:null, max:null });
    });
  }
  return cols;
}

// ═══════════════════════════════════════════════════════════════════
//  DESCRIPTIVE STATISTICS
// ═══════════════════════════════════════════════════════════════════

function descStats(arr) {
  const v = arr.filter(x => isFinite(x));
  if (!v.length) return { n:0, mean:NaN, sd:NaN, sem:NaN, median:NaN, q1:NaN, q3:NaN, min:NaN, max:NaN };
  const n    = v.length;
  const mean = v.reduce((a,b)=>a+b,0) / n;
  const sd   = n > 1 ? Math.sqrt(v.reduce((s,x)=>s+(x-mean)**2,0)/(n-1)) : 0;
  const sem  = sd / Math.sqrt(n);
  const sorted = [...v].sort((a,b)=>a-b);
  const median = _quantile(sorted, 0.5);
  const q1     = _quantile(sorted, 0.25);
  const q3     = _quantile(sorted, 0.75);
  return { n, mean, sd, sem, median, q1, q3, min:sorted[0], max:sorted[n-1] };
}

function _quantile(sorted, p) {
  const idx = p * (sorted.length - 1);
  const lo  = Math.floor(idx), hi = Math.ceil(idx);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (idx - lo);
}

// ═══════════════════════════════════════════════════════════════════
//  NORMALITY — SHAPIRO-WILK  (Royston 1995, n ≤ 5000)
// ═══════════════════════════════════════════════════════════════════

function shapiroWilk(x) {
  const v = [...x].filter(d=>isFinite(d)).sort((a,b)=>a-b);
  const n = v.length;
  if (n < 3)  return { W: NaN, p: NaN, note: "n < 3" };
  if (n > 5000) return { W: NaN, p: NaN, note: "n > 5000" };

  // Mean
  const mean = v.reduce((a,b)=>a+b,0)/n;

  // SS
  const ss = v.reduce((s,xi)=>s+(xi-mean)**2,0);

  // Coefficients via approximation (Royston)
  const m   = v.map((_,i) => jStat.normal.inv((i+1-0.375)/(n+0.25),0,1));
  const mss = m.reduce((s,mi)=>s+mi*mi,0);
  const c   = m.map(mi => mi/Math.sqrt(mss));

  // a coefficients (use last half)
  const half = Math.floor(n/2);
  let W_num  = 0;
  for (let i = 0; i < half; i++) {
    W_num += c[n-1-i] * (v[n-1-i] - v[i]);
  }
  const W = (W_num * W_num) / ss;

  // p-value via log(1-W) normal approximation (Royston 1995)
  const mu    = _swMu(n);
  const sigma = _swSigma(n);
  const z     = (Math.log(1 - W) - mu) / sigma;
  const p     = 1 - jStat.normal.cdf(z, 0, 1);

  return { W: +W.toFixed(4), p: +p.toFixed(4) };
}

function _swMu(n) {
  // Polynomial approximation of E[log(1-W)] from Royston 1995
  const ln = Math.log(n);
  return 0.0038915*(ln**3) - 0.083751*(ln**2) - 0.31082*ln - 1.5861;
}
function _swSigma(n) {
  const ln = Math.log(n);
  return Math.exp(0.0030302*(ln**2) - 0.082676*ln - 0.4803);
}

// ═══════════════════════════════════════════════════════════════════
//  TWO-SAMPLE TESTS
// ═══════════════════════════════════════════════════════════════════

// Welch's t-test (unpaired, unequal variance)
function tTestWelch(a, b, tails=2) {
  const da = descStats(a), db = descStats(b);
  if (da.n < 2 || db.n < 2) return null;
  const se  = Math.sqrt(da.sd**2/da.n + db.sd**2/db.n);
  if (se === 0) return { t:0, df:0, p:1, se };
  const t   = (da.mean - db.mean) / se;
  // Welch-Satterthwaite df
  const v1  = da.sd**2/da.n, v2 = db.sd**2/db.n;
  const df  = (v1+v2)**2 / (v1**2/(da.n-1) + v2**2/(db.n-1));
  const p1  = jStat.studentt.cdf(-Math.abs(t), df);
  const p   = tails === 1 ? p1 : 2*p1;
  const d   = cohensD(a, b);
  return { t:+t.toFixed(4), df:+df.toFixed(2), p:+p.toFixed(6), se:+se.toFixed(4), cohensD:d };
}

// Paired t-test
function tTestPaired(a, b, tails=2) {
  const diffs = a.map((ai,i) => ai - b[i]).filter(d=>isFinite(d));
  const n     = diffs.length;
  if (n < 2) return null;
  const ds    = descStats(diffs);
  const se    = ds.sd / Math.sqrt(n);
  if (se === 0) return { t:0, df:0, p:1 };
  const t     = ds.mean / se;
  const df    = n - 1;
  const p1    = jStat.studentt.cdf(-Math.abs(t), df);
  const p     = tails === 1 ? p1 : 2*p1;
  return { t:+t.toFixed(4), df, p:+p.toFixed(6), meanDiff:+ds.mean.toFixed(4), se:+se.toFixed(4) };
}

// Mann-Whitney U
function mannWhitneyU(a, b, tails=2) {
  const na = a.length, nb = b.length;
  if (!na || !nb) return null;
  let U1 = 0;
  a.forEach(ai => b.forEach(bi => {
    if (ai > bi) U1++;
    else if (ai === bi) U1 += 0.5;
  }));
  const U2 = na*nb - U1;
  const U  = Math.min(U1, U2);
  // Normal approximation with tie correction
  const mu_U = na*nb/2;
  // Tie correction
  const combined = [...a,...b].sort((x,y)=>x-y);
  let tieSum = 0;
  let i = 0;
  while (i < combined.length) {
    let j = i;
    while (j < combined.length && combined[j]===combined[i]) j++;
    const t = j - i;
    if (t > 1) tieSum += t**3 - t;
    i = j;
  }
  const N = na + nb;
  const sigma_U = Math.sqrt((na*nb/12)*((N+1) - tieSum/(N*(N-1))));
  if (sigma_U === 0) return { U, p:1, r:0 };
  const z  = (U - mu_U) / sigma_U;
  const p1 = jStat.normal.cdf(z, 0, 1);  // U is min so z is negative
  const p  = tails === 1 ? p1 : 2*Math.min(p1, 1-p1);
  const r  = Math.abs(z) / Math.sqrt(N);   // rank-biserial effect size
  return { U:+U.toFixed(2), z:+z.toFixed(4), p:+Math.max(0,p).toFixed(6), r:+r.toFixed(4) };
}

// ═══════════════════════════════════════════════════════════════════
//  MULTI-GROUP TESTS
// ═══════════════════════════════════════════════════════════════════

// One-way ANOVA (returns F, df, p, eta-squared)
function oneWayAnova(groups) {
  // groups: { label: string, values: number[] }[]
  const k    = groups.length;
  const ns   = groups.map(g=>g.values.length);
  const N    = ns.reduce((a,b)=>a+b,0);
  if (N < k + 1 || k < 2) return null;

  const grandMean = groups.flatMap(g=>g.values).reduce((a,b)=>a+b,0) / N;
  const SSb = groups.reduce((s,g)=>s + g.values.length*(descStats(g.values).mean - grandMean)**2, 0);
  const SSw = groups.reduce((s,g)=>{ const m=descStats(g.values).mean; return s+g.values.reduce((ss,v)=>ss+(v-m)**2,0);},0);
  const dfb = k - 1, dfw = N - k;
  if (dfw <= 0 || SSw === 0) return null;
  const F   = (SSb/dfb) / (SSw/dfw);
  const p   = 1 - jStat.centralF.cdf(F, dfb, dfw);
  const eta2 = SSb / (SSb + SSw);
  return { F:+F.toFixed(4), dfb, dfw, p:+p.toFixed(6), eta2:+eta2.toFixed(4), SSb, SSw };
}

// Kruskal-Wallis H test
function kruskalWallis(groups) {
  const k = groups.length;
  if (k < 2) return null;
  const allVals = groups.flatMap(g=>g.values);
  const N = allVals.length;
  // Rank all values
  const sorted = [...allVals].sort((a,b)=>a-b);
  const rankMap = new Map();
  let i = 0;
  while (i < sorted.length) {
    let j = i;
    while (j < sorted.length && sorted[j]===sorted[i]) j++;
    const avgRank = (i + j + 1) / 2;   // 1-based
    for (let k2=i; k2<j; k2++) {
      const v = sorted[k2];
      if (!rankMap.has(v)) rankMap.set(v,[]);
      rankMap.get(v).push(avgRank);
    }
    i = j;
  }
  const getRank = v => { const rs=rankMap.get(v)||[]; return rs.reduce((a,b)=>a+b,0)/rs.length; };

  // Tie correction
  let C = 0;
  rankMap.forEach((_,v) => {
    const t = rankMap.get(v).length;
    C += t**3 - t;
  });
  const tieCorr = 1 - C/(N**3-N);

  const H_num = groups.reduce((s,g)=>{
    const Rj = g.values.reduce((a,v)=>a+getRank(v),0);
    return s + Rj**2/g.values.length;
  },0);
  const H_raw = (12/(N*(N+1))) * H_num - 3*(N+1);
  const H = tieCorr > 0 ? H_raw / tieCorr : H_raw;
  const df = k - 1;
  const p  = 1 - jStat.chisquare.cdf(H, df);
  return { H:+H.toFixed(4), df, p:+p.toFixed(6) };
}

// ═══════════════════════════════════════════════════════════════════
//  POST-HOC TESTS
// ═══════════════════════════════════════════════════════════════════

// Returns array of { groupA, groupB, p_raw, p_adj, significant }
function postHocTukey(groups) {
  const k  = groups.length;
  const ns = groups.map(g=>g.values.length);
  const N  = ns.reduce((a,b)=>a+b,0);
  const ds = groups.map(g=>descStats(g.values));
  const MSw = groups.reduce((s,g,i)=>{
    const m=ds[i].mean; return s+g.values.reduce((ss,v)=>ss+(v-m)**2,0);
  },0) / (N-k);

  const pairs = [];
  for (let i=0;i<k;i++) for (let j=i+1;j<k;j++) {
    const diff  = Math.abs(ds[i].mean - ds[j].mean);
    const se    = Math.sqrt(MSw/2*(1/ns[i]+1/ns[j]));
    if (se === 0) { pairs.push({groupA:groups[i].label, groupB:groups[j].label, p_raw:1, p_adj:1}); continue; }
    const q    = diff / se;
    // Studentized range p-value via jStat (approximation: treat q/√2 as t)
    const t    = q / Math.sqrt(2);
    const df   = N - k;
    const p_raw = 2*jStat.studentt.cdf(-Math.abs(t), df);
    // Tukey adjustment via Bonferroni approximation on k(k-1)/2 comparisons
    const m    = k*(k-1)/2;
    const p_adj = Math.min(1, p_raw * m);
    pairs.push({ groupA:groups[i].label, groupB:groups[j].label, p_raw:+p_raw.toFixed(6), p_adj:+p_adj.toFixed(6) });
  }
  return pairs;
}

// Dunn's test (post-hoc for Kruskal-Wallis)
function postHocDunn(groups) {
  const allVals = groups.flatMap(g=>g.values);
  const N = allVals.length;
  const sorted = [...allVals].sort((a,b)=>a-b);

  // Assign average ranks
  const rankOf = new Map();
  let i2=0;
  while (i2<sorted.length) {
    let j=i2;
    while(j<sorted.length && sorted[j]===sorted[i2]) j++;
    const avg=(i2+j+1)/2;
    for(let k=i2;k<j;k++) { const v=sorted[k]; if(!rankOf.has(v)) rankOf.set(v,avg); }
    i2=j;
  }

  // Tie correction term
  const tieGroups = new Map();
  sorted.forEach(v=>{tieGroups.set(v,(tieGroups.get(v)||0)+1);});
  let T=0; tieGroups.forEach(t=>{ if(t>1) T+=t**3-t; });
  const tieVar = (N*(N+1)/12) - T/(12*(N-1));

  const groupRankMeans = groups.map(g=>{
    const Rj = g.values.reduce((s,v)=>s+rankOf.get(v),0);
    return Rj / g.values.length;
  });

  const pairs=[];
  const k=groups.length;
  for(let a=0;a<k;a++) for(let b=a+1;b<k;b++){
    const se = Math.sqrt(tieVar*(1/groups[a].values.length+1/groups[b].values.length));
    if(se===0){pairs.push({groupA:groups[a].label,groupB:groups[b].label,p_raw:1,p_adj:1});continue;}
    const z = Math.abs(groupRankMeans[a]-groupRankMeans[b])/se;
    const p_raw = 2*(1-jStat.normal.cdf(z,0,1));
    pairs.push({groupA:groups[a].label,groupB:groups[b].label,p_raw:+p_raw.toFixed(6),p_adj:p_raw});
  }
  return pairs;
}

// Apply multiple-comparison corrections to an array of pair objects
function applyCorrection(pairs, method) {
  if (method === "none") return pairs.map(p=>({...p, p_adj:p.p_raw}));
  if (method === "bonferroni") {
    const m = pairs.length;
    return pairs.map(p=>({...p, p_adj:+Math.min(1,p.p_raw*m).toFixed(6)}));
  }
  if (method === "bh") {
    const m = pairs.length;
    const sorted = [...pairs].sort((a,b)=>a.p_raw-b.p_raw);
    const adj = sorted.map((p,i)=>({...p, p_adj:+Math.min(1,p.p_raw*m/(i+1)).toFixed(6)}));
    // Enforce monotonicity
    for(let i=adj.length-2;i>=0;i--) adj[i].p_adj=Math.min(adj[i].p_adj,adj[i+1].p_adj);
    // Re-map back to original order
    const adjMap=new Map(adj.map(p=>[p.groupA+"__"+p.groupB,p.p_adj]));
    return pairs.map(p=>({...p,p_adj:adjMap.get(p.groupA+"__"+p.groupB)??p.p_adj}));
  }
  return pairs;
}

// ═══════════════════════════════════════════════════════════════════
//  CORRELATION & REGRESSION
// ═══════════════════════════════════════════════════════════════════

function pearsonCorr(x, y) {
  const n = Math.min(x.length, y.length);
  if (n < 3) return null;
  const mx = x.slice(0,n).reduce((a,b)=>a+b,0)/n;
  const my = y.slice(0,n).reduce((a,b)=>a+b,0)/n;
  let num=0,dx2=0,dy2=0;
  for(let i=0;i<n;i++){
    const ex=x[i]-mx, ey=y[i]-my;
    num+=ex*ey; dx2+=ex**2; dy2+=ey**2;
  }
  const r  = num/Math.sqrt(dx2*dy2);
  const t  = r*Math.sqrt((n-2)/(1-r**2));
  const p  = 2*jStat.studentt.cdf(-Math.abs(t),n-2);
  return { r:+r.toFixed(4), r2:+(r**2).toFixed(4), t:+t.toFixed(4), df:n-2, p:+p.toFixed(6), method:"Pearson" };
}

function spearmanCorr(x, y) {
  const n = Math.min(x.length,y.length);
  if (n < 3) return null;
  const rank = arr => {
    const sorted=[...arr].map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);
    const ranks=new Array(n);
    let i2=0;
    while(i2<n){
      let j=i2; while(j<n&&sorted[j].v===sorted[i2].v) j++;
      const avg=(i2+j+1)/2;
      for(let k=i2;k<j;k++) ranks[sorted[k].i]=avg;
      i2=j;
    }
    return ranks;
  };
  const rx=rank(x.slice(0,n)), ry=rank(y.slice(0,n));
  return pearsonCorr(rx,ry) ? {...pearsonCorr(rx,ry), method:"Spearman"} : null;
}

function linearRegression(x, y) {
  const n = Math.min(x.length,y.length);
  if (n < 2) return null;
  const mx=x.slice(0,n).reduce((a,b)=>a+b,0)/n;
  const my=y.slice(0,n).reduce((a,b)=>a+b,0)/n;
  let Sxy=0,Sxx=0;
  for(let i=0;i<n;i++){Sxy+=(x[i]-mx)*(y[i]-my);Sxx+=(x[i]-mx)**2;}
  if(Sxx===0) return null;
  const slope=Sxy/Sxx, intercept=my-slope*mx;
  const yhat=x.slice(0,n).map(xi=>slope*xi+intercept);
  const SSres=y.slice(0,n).reduce((s,yi,i)=>s+(yi-yhat[i])**2,0);
  const SStot=y.slice(0,n).reduce((s,yi)=>s+(yi-my)**2,0);
  const r2=SStot>0?1-SSres/SStot:1;
  return { slope:+slope.toFixed(6), intercept:+intercept.toFixed(6), r2:+r2.toFixed(4) };
}

// ═══════════════════════════════════════════════════════════════════
//  EFFECT SIZES
// ═══════════════════════════════════════════════════════════════════

function cohensD(a, b) {
  const da=descStats(a), db=descStats(b);
  const pooledSD = Math.sqrt(((da.n-1)*da.sd**2+(db.n-1)*db.sd**2)/(da.n+db.n-2));
  if(pooledSD===0) return 0;
  return +((da.mean-db.mean)/pooledSD).toFixed(4);
}

// ═══════════════════════════════════════════════════════════════════
//  SIGNIFICANCE LABEL
// ═══════════════════════════════════════════════════════════════════

function sigStars(p) {
  if(p===null||p===undefined||isNaN(p)) return "ns";
  if(p < 0.0001) return "****";
  if(p < 0.001)  return "***";
  if(p < 0.01)   return "**";
  if(p < 0.05)   return "*";
  return "ns";
}

function fmtP(p, mode="exact") {
  if(p===null||p===undefined||isNaN(p)) return "ns";
  if(mode==="stars") return sigStars(p);
  if(p < 0.0001) return "p < 0.0001";
  return "p = " + p.toFixed(4).replace(/0+$/,"").replace(/\.$/,"");
}
