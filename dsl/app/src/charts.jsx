import React from "react";

const scale = (d0, d1, r0, r1) => (v) => r0 + ((v - d0) / (d1 - d0)) * (r1 - r0);
export const ACC = [0.17, 0.47]; // the shared accuracy axis for every decode chart

export function ticks(lo, hi, n = 6) {
  const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / n)));
  const err = (hi - lo) / n / step;
  const s = step * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1);
  const out = [];
  for (let v = Math.ceil(lo / s) * s; v <= hi + 1e-12; v += s) out.push(+v.toFixed(10));
  return out;
}

// Student t 97.5% quantile for small df (two-sided 95% interval).
const T975 = [0, 12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.36, 2.31, 2.26, 2.23, 2.20, 2.18, 2.16, 2.14, 2.13, 2.12, 2.11, 2.10, 2.09, 2.09];
export function meanCI(v) {
  const n = v.length, m = v.reduce((a, b) => a + b, 0) / n;
  const sd = Math.sqrt(v.reduce((a, b) => a + (b - m) ** 2, 0) / (n - 1));
  const h = (T975[n - 1] || 1.96) * sd / Math.sqrt(n);
  return { n, mean: m, lo: m - h, hi: m + h };
}

function Axis({ x, lo, hi, y, top }) {
  return ticks(lo, hi).map((t) => (
    <g key={t}>
      <line x1={x(t)} x2={x(t)} y1={top} y2={y} stroke="var(--line-2)" />
      <text x={x(t)} y={y + 16} textAnchor="middle" className="tnum">{t.toFixed(2)}</text>
    </g>
  ));
}

// One subject: its searchlight centres against its own pooled null, as
// probability densities on one square-root scale. The FDR threshold is drawn
// with the null's tail beyond it shaded, and the surviving centres as a rug.
export function SubjectNull({ subject, edges, chance, resolution, permutations, domain = ACC, width = 1060 }) {
  const [lo, hi] = domain;
  const at = (f) => lo + f * (hi - lo);
  const m = { l: 120, r: 24 };
  const x = scale(lo, hi, m.l, width - m.r);
  const H = 150, base = 120;
  const dens = (c) => { const s = c.reduce((a, b) => a + b, 0); return c.map((v) => v / s); };
  const nd = dens(subject.null_hist), od = dens(subject.centres_hist);
  const ymax = Math.sqrt(Math.max(...nd, ...od));
  const y = (d) => base - (Math.sqrt(d) / ymax) * (base - 14);
  const inDomain = (i) => edges[i] >= lo && edges[i + 1] <= hi;
  const step = (d) => d.map((v, i) => (inDomain(i) ? [i, v] : null)).filter(Boolean)
    .map(([i, v], k) => `${k ? "L" : "M"}${x(edges[i]).toFixed(1)},${y(v).toFixed(1)}L${x(edges[i + 1]).toFixed(1)},${y(v).toFixed(1)}`).join("");
  const thr = subject.fdr?.threshold;
  const tail = thr ? nd.map((v, i) => (edges[i] >= thr ? `M${x(edges[i])},${base}L${x(edges[i])},${y(v)}L${x(Math.min(hi, edges[i + 1]))},${y(v)}L${x(Math.min(hi, edges[i + 1]))},${base}Z` : "")).join("") : "";
  const rug = thr ? od.map((v, i) => (edges[i] >= thr && v > 0 && inDomain(i) ? x(edges[i] + (edges[1] - edges[0]) / 2) : null)).filter((v) => v != null) : [];
  return (
    <svg className="chart" width={width} height={H + 26} role="img" aria-label={`sub-${subject.id}: centres against its own null`}>
      <text x={0} y={22} style={{ fontWeight: 600, fill: "var(--ink)" }}>sub-{subject.id}</text>
      <text x={0} y={38}>its centres vs</text>
      <text x={0} y={52}>its own null</text>
      <text x={0} y={70} style={{ fontSize: 10 }}>√ density</text>
      <Axis x={x} lo={lo} hi={hi} y={base} top={8} />
      <path d={step(nd) + `L${x(hi)},${base}L${x(edges.find((e) => e >= lo))},${base}Z`} fill="var(--ink-3)" opacity={0.25} />
      <path d={tail} fill="var(--warm)" opacity={0.35} />
      <path d={step(od)} fill="none" stroke="var(--ink)" strokeWidth={1.4} />
      <line x1={x(chance)} x2={x(chance)} y1={6} y2={base} stroke="var(--ink)" strokeDasharray="2 3" />
      <text x={x(chance) - 4} y={16} textAnchor="end" style={{ fill: "var(--ink)" }}>chance</text>
      <text x={x(at(0.12))} y={y(Math.max(...nd)) + 30} textAnchor="end" style={{ fill: "var(--ink-2)" }}>null</text>
      <text x={x(at(0.4))} y={y(od[Math.round((at(0.4) - edges[0]) / (edges[1] - edges[0]))] || 0) - 6} style={{ fill: "var(--ink)" }}>centres</text>
      {thr ? <>
        <line x1={x(thr)} x2={x(thr)} y1={6} y2={base + 8} stroke="var(--warm)" strokeWidth={2} />
        <text x={x(thr) + 6} y={18} style={{ fill: "var(--warm)", fontWeight: 600 }}>FDR q &lt; .05 at ≥ {thr.toFixed(3)}{permutations < 1000 ? " · provisional" : ""}</text>
        <text x={x(thr) + 6} y={33} style={{ fill: "var(--warm)" }}>{subject.fdr.surviving.toLocaleString()} centres survive</text>
        {rug.map((rx, i) => <line key={i} x1={rx} x2={rx} y1={base - 8} y2={base - 1} stroke="var(--warm)" strokeWidth={1.5} />)}
      </> : <text x={x(at(0.63))} y={24} style={{ fill: "var(--ink-2)", fontWeight: 600 }}>no centre survives FDR q &lt; .05</text>}
      <text x={width - m.r} y={H + 20} textAnchor="end" style={{ fill: "var(--ink-3)" }}>
        y: square root of density · p resolution ≈ {resolution} with {permutations} permutations{permutations < 1000 ? "; counts are provisional below 1,000" : ""}</text>
    </svg>
  );
}

// Subjects on the same axis: accuracy in a leave-one-subject-out region,
// the mean with a 95% t-interval, extremes labelled, stale shown hollow.
export function SubjectStrip({ subjects, chance, include, hover, setHover, select, selected, domain = ACC, width = 1060 }) {
  const [lo, hi] = domain;
  const m = { l: 120, r: 24 };
  const x = scale(lo, hi, m.l, width - m.r);
  const subs = subjects.filter((s) => typeof s.loso_roi === "number");
  const used = subs.filter((s) => s.status === "done" || include);
  const ci = meanCI(used.map((s) => s.loso_roi));
  const sorted = [...used].sort((a, b) => a.loso_roi - b.loso_roi);
  const extremes = new Set([sorted[0]?.id, sorted[sorted.length - 1]?.id]);
  // Centred 1-D beeswarm: rows alternate above and below one line.
  const placed = [];
  const order = [0, -1, 1, -2, 2, -3, 3];
  const pts = [...subs].sort((a, b) => a.loso_roi - b.loso_roi).map((s) => {
    const px = x(s.loso_roi);
    const row = order.find((r) => !placed.some((p) => p.row === r && Math.abs(p.px - px) < 11)) ?? 4;
    placed.push({ px, row });
    return { s, px, py: 44 + row * 10 };
  });
  const H = 112;
  return (
    <svg className="chart" width={width} height={H + 20} role="img" aria-label="Each subject's accuracy in a leave-one-subject-out region">
      <text x={0} y={22} style={{ fontWeight: 600, fill: "var(--ink)" }}>Subjects</text>
      <text x={0} y={38}>accuracy in the</text>
      <text x={0} y={52}>LOSO region</text>
      <Axis x={x} lo={lo} hi={hi} y={H} top={8} />
      <line x1={x(chance)} x2={x(chance)} y1={6} y2={H} stroke="var(--ink)" strokeDasharray="2 3" />
      <rect x={x(ci.lo)} y={H - 26} width={x(ci.hi) - x(ci.lo)} height={8} fill="var(--ink)" opacity={0.14} />
      <line x1={x(ci.mean)} x2={x(ci.mean)} y1={H - 30} y2={H - 14} stroke="var(--ink)" strokeWidth={2} />
      <text x={x(ci.hi) + 6} y={H - 18} style={{ fill: "var(--ink)" }}>mean {ci.mean.toFixed(3)}, 95% CI {ci.lo.toFixed(3)}–{ci.hi.toFixed(3)} · n = {ci.n}</text>
      {pts.map(({ s, px, py }) => {
        const stale = s.status !== "done";
        const on = hover === s.id;
        const lab = on || extremes.has(s.id) || stale;
        return (
          <g key={s.id} onMouseEnter={() => setHover(s.id)} onMouseLeave={() => setHover(null)} onClick={() => select?.(s.id)} style={{ cursor: "pointer" }}>
            {selected === s.id && <circle cx={px} cy={py} r={8.5} fill="none" stroke="var(--ink)" strokeWidth={2} />}
            <circle cx={px} cy={py} r={on ? 6 : 4.5} fill={stale ? "var(--surface)" : "var(--ink)"} stroke={stale ? "var(--amber)" : "var(--surface)"} strokeWidth={stale ? 1.8 : 1} />
            {stale && <line x1={px} x2={px + 18} y1={py + 5} y2={py + 22} stroke="var(--amber)" />}
            {lab && <text x={stale ? px + 21 : px} y={stale ? py + 26 : py - 10} textAnchor={stale ? "start" : "middle"} style={{ fill: stale ? "var(--amber)" : "var(--ink)", fontWeight: on ? 600 : 400 }}>
              {s.id}{stale ? (include ? " stale, in mean" : " stale, not in mean") : on ? ` ${s.loso_roi.toFixed(3)}` : ""}</text>}
          </g>
        );
      })}
    </svg>
  );
}

// Keyed colour ramp on the accuracy axis.
export function AccuracyKey({ chance, domain = ACC, width = 420, label = "own FDR threshold", mark }) {
  const [lo, hi] = domain;
  const x = scale(lo, hi, 0, width);
  const t = mark ?? chance + 0.45 * (hi - chance);
  return (
    <svg className="chart" width={width + 20} height={34} role="img" aria-label="Colour key">
      <defs>
        <linearGradient id="faint" x1="0" x2="1"><stop offset="0" stopColor="#d9d6d0" /><stop offset="1" stopColor="#e8b9ac" /></linearGradient>
        <linearGradient id="hot" x1="0" x2="1"><stop offset="0" stopColor="#7a1d0e" /><stop offset=".35" stopColor="#d4492c" /><stop offset=".7" stopColor="#f3b13b" /><stop offset="1" stopColor="#fff1c2" /></linearGradient>
      </defs>
      <rect x={x(chance)} y={4} width={x(t) - x(chance)} height={8} fill="url(#faint)" />
      <rect x={x(t)} y={4} width={x(hi) - x(t)} height={8} fill="url(#hot)" />
      <text x={x(chance)} y={28} textAnchor="middle">chance</text>
      <text x={x(t)} y={28} textAnchor="middle">{label}</text>
      <text x={x(hi)} y={28} textAnchor="end">{hi.toFixed(2)}</text>
    </svg>
  );
}

// Regions, paired by hemisphere: one row per region pair, LH dots filled and
// RH dots hollow on a slight vertical offset, a mean tick and 95% interval per
// hemisphere. Hover links a subject across rows.
export function RegionPairs({ regions, subjects, chance, order, hover, setHover, width = 1060 }) {
  const half = regions.length / 2;
  const stat = (r) => { const v = r.values.map((x) => (typeof x === "number" ? x : null)); return { ...r, v, ...meanCI(v.filter((x) => x != null)) }; };
  let pairs = Array.from({ length: half }, (_, i) => ({ L: stat(regions[i]), R: stat(regions[i + half]), name: regions[i].name.replace(/^LH /, "") }));
  if (order === "mean") pairs.sort((a, b) => Math.max(b.L.mean, b.R.mean) - Math.max(a.L.mean, a.R.mean));
  const all = pairs.flatMap((p) => [...p.L.v, ...p.R.v].filter((x) => x != null));
  const lo = Math.min(chance, ...all) - 0.01, hi = Math.ceil((Math.max(...all) + 0.005) / 0.05) * 0.05;
  const rowH = 30, m = { l: 210, r: 190, t: 26, b: 26 };
  const H = m.t + pairs.length * rowH + m.b;
  const x = scale(lo, hi, m.l, width - m.r);
  const jitter = (j) => ((j * 7919) % 5 - 2) * 1.1;
  const Hemi = ({ h, dy, filled }) => (
    <g transform={`translate(0, ${dy})`}>
      {h.v.map((v, j) => v == null ? null : (
        <circle key={j} cx={x(v)} cy={jitter(j)} r={hover === subjects[j] ? 4.5 : 2.4}
          fill={filled ? "var(--ink)" : "var(--surface)"} stroke="var(--ink)" strokeWidth={filled ? 0 : 1}
          opacity={hover == null ? 0.45 : hover === subjects[j] ? 1 : 0.12} onMouseEnter={() => setHover(subjects[j])} />
      ))}
      <line x1={x(h.lo)} x2={x(h.hi)} y1={0} y2={0} stroke="var(--ink)" strokeWidth={2} />
      <line x1={x(h.mean)} x2={x(h.mean)} y1={-6} y2={6} stroke="var(--ink)" strokeWidth={2.5} />
    </g>
  );
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label="Accuracy by region pair, every subject" onMouseLeave={() => setHover(null)}>
      {ticks(lo, hi).map((t) => (
        <g key={t}>
          <line x1={x(t)} x2={x(t)} y1={m.t - 6} y2={H - m.b} stroke="var(--line-2)" />
          <text x={x(t)} y={H - 8} textAnchor="middle">{t.toFixed(2)}</text>
        </g>
      ))}
      <line x1={x(chance)} x2={x(chance)} y1={m.t - 10} y2={H - m.b} stroke="var(--ink)" strokeDasharray="2 3" />
      <text x={x(chance) + 4} y={m.t - 12} style={{ fill: "var(--ink)" }}>chance {chance}</text>
      <text x={width - m.r + 12} y={m.t - 12}>LH ● · RH ○ mean (95% CI)</text>
      {hover && <text x={m.l} y={m.t - 12} style={{ fill: "var(--ink)", fontWeight: 600 }}>sub-{hover} in every row</text>}
      {pairs.map((p, i) => (
        <g key={p.name} transform={`translate(0, ${m.t + i * rowH})`}>
          <text x={m.l - 12} y={rowH / 2 + 4} textAnchor="end">{p.name}</text>
          <Hemi h={p.L} dy={rowH / 2 - 6} filled />
          <Hemi h={p.R} dy={rowH / 2 + 6} />
          <text x={width - m.r + 12} y={rowH / 2 - 2} className="tnum">● {p.L.mean.toFixed(3)} <tspan style={{ fill: "var(--ink-3)" }}>{p.L.lo.toFixed(2)}–{p.L.hi.toFixed(2)}</tspan></text>
          <text x={width - m.r + 12} y={rowH / 2 + 12} className="tnum">○ {p.R.mean.toFixed(3)} <tspan style={{ fill: "var(--ink-3)" }}>{p.R.lo.toFixed(2)}–{p.R.hi.toFixed(2)}</tspan></text>
        </g>
      ))}
    </svg>
  );
}

// Binned key for subjects-significant counts; the same bins colour the map.
export function CountKey({ bins }) {
  return (
    <svg className="chart" width={bins.length * 46} height={34} role="img" aria-label="Count key">
      {bins.map((b, i) => (
        <g key={b.label}>
          <rect x={i * 46} y={4} width={42} height={10} fill={`rgb(${b.rgb.join(",")})`} />
          <text x={i * 46 + 21} y={28} textAnchor="middle">{b.label}</text>
        </g>
      ))}
    </svg>
  );
}
