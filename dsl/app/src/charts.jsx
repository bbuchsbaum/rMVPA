import React from "react";

const scale = (d0, d1, r0, r1) => (v) => r0 + ((v - d0) / (d1 - d0)) * (r1 - r0);

export function ticks(lo, hi, n = 6) {
  const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / n)));
  const err = (hi - lo) / n / step;
  const s = step * (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1);
  const out = [];
  for (let v = Math.ceil(lo / s) * s; v <= hi + 1e-12; v += s) out.push(+v.toFixed(10));
  return out;
}

// Two panels on one accuracy axis:
//  top    — every searchlight centre (warm outline) against the pooled
//           permutation null (grey), each scaled to its own peak, with the
//           FDR threshold and chance drawn and labelled;
//  bottom — one dot per subject: its accuracy in a leave-one-subject-out
//           region. Stale subjects are hollow amber; hover links to the tiles.
export function SharedAxis({ res, include, hover, setHover, width = 1060 }) {
  const lo = 0.17, hi = 0.45;
  const m = { l: 120, r: 24 };
  const x = scale(lo, hi, m.l, width - m.r);
  const { edges } = res.null.hist;
  const nullC = res.null.hist.counts, obsC = res.centres.hist.counts;
  const top = 120, gap = 34, bottom = 96, H = top + gap + bottom + 26;
  const yN = scale(0, Math.max(...nullC), top, 12), yO = scale(0, Math.max(...obsC), top, 12);
  const step = (counts, y) => counts.map((c, i) => `${i ? "L" : "M"}${x(edges[i]).toFixed(1)},${y(c).toFixed(1)}L${x(edges[i + 1]).toFixed(1)},${y(c).toFixed(1)}`).join("");
  const fdr = res.fdr.threshold;
  const subs = res.subjects.filter((s) => s.loso_roi != null);
  const used = subs.filter((s) => s.status === "done" || include);
  const mean = used.reduce((a, s) => a + s.loso_roi, 0) / used.length;
  const y0 = top + gap;
  // Beeswarm-free: dots on one row with a deterministic vertical offset only when they would overlap.
  const placed = [];
  const pts = subs.map((s) => {
    const px = x(s.loso_roi);
    let row = 0;
    while (placed.some((p) => p.row === row && Math.abs(p.px - px) < 11)) row++;
    placed.push({ px, row });
    return { s, px, py: y0 + 22 + row * 12 };
  });
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label="Searchlight centres and subjects against the null">
      <text x={0} y={22} style={{ fontWeight: 600, fill: "var(--ink)" }}>All centres</text>
      <text x={0} y={38}>{res.centres.n.toLocaleString()} searchlights</text>
      <text x={0} y={y0 + 16} style={{ fontWeight: 600, fill: "var(--ink)" }}>Subjects</text>
      <text x={0} y={y0 + 32}>LOSO region</text>
      {ticks(lo, hi).map((t) => (
        <g key={t}>
          <line x1={x(t)} x2={x(t)} y1={10} y2={H - 26} stroke="var(--line-2)" />
          <text x={x(t)} y={H - 8} textAnchor="middle">{t.toFixed(2)}</text>
        </g>
      ))}
      <path d={step(nullC, yN) + `L${x(edges[edges.length - 1])},${top}L${x(edges[0])},${top}Z`} fill="var(--ink-3)" opacity={0.28} />
      <path d={step(obsC, yO)} fill="none" stroke="var(--warm)" strokeWidth={1.5} />
      <text x={x(0.272)} y={yN(Math.max(...nullC)) + 34} style={{ fill: "var(--ink-2)", fontWeight: 600 }}>pooled null</text>
      <text x={x(0.272)} y={yN(Math.max(...nullC)) + 48}>{res.null.permutations} perms × 10% of centres</text>
      <text x={x(0.31)} y={yO(obsC[Math.floor((0.31 - edges[0]) / (edges[1] - edges[0]))] || 0) - 6} style={{ fill: "var(--warm)" }}>observed centres</text>
      <line x1={x(res.chance)} x2={x(res.chance)} y1={6} y2={H - 26} stroke="var(--ink)" strokeDasharray="2 3" />
      <text x={x(res.chance) + 4} y={H - 30} style={{ fill: "var(--ink)" }}>chance {res.chance}</text>
      {fdr && <>
        <line x1={x(fdr)} x2={x(fdr)} y1={6} y2={top} stroke="var(--warm)" strokeWidth={2} />
        <text x={x(fdr) + 6} y={18} style={{ fill: "var(--warm)", fontWeight: 600 }}>FDR q &lt; .05 → ≥ {fdr.toFixed(3)}</text>
        <text x={x(fdr) + 6} y={33} style={{ fill: "var(--warm)" }}>{res.fdr.surviving} of {res.centres.n.toLocaleString()} centres</text>
      </>}
      <line x1={m.l} x2={width - m.r} y1={y0} y2={y0} stroke="var(--line)" />
      {pts.map(({ s, px, py }) => {
        const stale = s.status !== "done";
        const on = hover === s.id;
        return (
          <g key={s.id} onMouseEnter={() => setHover(s.id)} onMouseLeave={() => setHover(null)} style={{ cursor: "default" }}>
            <circle cx={px} cy={py} r={on ? 6 : 4.5}
              fill={stale ? "var(--surface)" : "var(--ink)"} stroke={stale ? "var(--amber)" : "var(--surface)"} strokeWidth={stale ? 1.8 : 1}
              opacity={stale && !include ? 0.9 : 1} />
            {on && <text x={px} y={py - 9} textAnchor="middle" style={{ fill: "var(--ink)", fontWeight: 600 }}>sub-{s.id} {s.loso_roi.toFixed(3)}{stale ? " (stale)" : ""}</text>}
          </g>
        );
      })}
      <line x1={x(mean)} x2={x(mean)} y1={y0 + 8} y2={y0 + bottom - 18} stroke="var(--ink)" strokeWidth={2} />
      <text x={x(mean) + 6} y={y0 + bottom - 20} style={{ fill: "var(--ink)" }}>mean {mean.toFixed(3)} · n = {used.length}</text>
    </svg>
  );
}

// Regions: one row per region, every subject as a dot, a mean tick, chance.
export function RegionDots({ regions, subjects, chance, order, width = 1060 }) {
  const rows = regions.map((r) => {
    const v = r.values.filter((x) => typeof x === "number" && x === x);
    return { ...r, v, mean: v.reduce((a, b) => a + b, 0) / v.length };
  });
  if (order === "mean") rows.sort((a, b) => b.mean - a.mean);
  const all = rows.flatMap((r) => r.v);
  const lo = Math.min(chance, ...all) - 0.01, hi = Math.max(...all) + 0.01;
  const rowH = 19, m = { l: 210, r: 70, t: 22, b: 26 };
  const H = m.t + rows.length * rowH + m.b;
  const x = scale(lo, hi, m.l, width - m.r);
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label="Accuracy by region, every subject">
      {ticks(lo, hi).map((t) => (
        <g key={t}>
          <line x1={x(t)} x2={x(t)} y1={m.t - 6} y2={H - m.b} stroke="var(--line-2)" />
          <text x={x(t)} y={H - 8} textAnchor="middle">{t.toFixed(2)}</text>
        </g>
      ))}
      <line x1={x(chance)} x2={x(chance)} y1={m.t - 10} y2={H - m.b} stroke="var(--ink)" strokeDasharray="2 3" />
      <text x={x(chance) + 4} y={m.t - 10} style={{ fill: "var(--ink)" }}>chance {chance}</text>
      <text x={width - m.r + 8} y={m.t - 10}>mean</text>
      {rows.map((r, i) => (
        <g key={r.index} transform={`translate(0, ${m.t + i * rowH})`}>
          {i % 2 === 0 && <rect x={0} y={0} width={width} height={rowH} fill="var(--sunk)" opacity={0.45} />}
          <text x={m.l - 10} y={13} textAnchor="end"><tspan style={{ fill: "var(--ink-3)" }}>{String(r.index).padStart(2, " ")}  </tspan>{r.name}</text>
          {r.v.map((v, j) => <circle key={j} cx={x(v)} cy={rowH / 2} r={2.6} fill="var(--ink)" opacity={0.35} />)}
          <line x1={x(r.mean)} x2={x(r.mean)} y1={3} y2={rowH - 3} stroke="var(--ink)" strokeWidth={2} />
          <text x={width - m.r + 8} y={13} className="tnum">{r.mean.toFixed(3)}</text>
        </g>
      ))}
    </svg>
  );
}
