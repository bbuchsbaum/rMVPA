// Whole-brain results: one value per subject per metric, on one axis each,
// with the mean and its 95% t-interval. Descriptive, like the other views.
import React, { useState } from "react";
import { meanCI, ticks } from "./charts.jsx";

const scale = (d0, d1, r0, r1) => (v) => r0 + ((v - d0) / (d1 - d0)) * (r1 - r0);

function MetricStrip({ metric, rows, chance, width = 1060 }) {
  const vals = rows.map((r) => r.v);
  const st = meanCI(vals);
  const lo = Math.min(chance, ...vals), hi = Math.max(chance, ...vals);
  const pad = (hi - lo || 1) * 0.08;
  const [d0, d1] = [lo - pad, hi + pad];
  const x = scale(d0, d1, 120, width - 230);
  // Dots that would overlap stack upward (a one-sided beeswarm), so four
  // subjects at 1.00 read as four, not one.
  const placed = [];
  const dots = [...rows].sort((p, q) => p.v - q.v).map((r) => {
    const cx = x(r.v);
    let level = 0;
    while (placed.some((d) => d.level === level && Math.abs(d.cx - cx) < 8)) level++;
    const d = { ...r, cx, level }; placed.push(d); return d;
  });
  const top = Math.max(0, ...dots.map((d) => d.level));
  const y = 30 + top * 8, H = y + 36;
  const m = { l: 120, r: 230 };
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label={`${metric}: every subject`}>
      <text x={0} y={y + 4} style={{ fontWeight: 600, fill: "var(--ink)" }}>{metric}</text>
      {ticks(d0, d1).map((t) => <g key={t}><line x1={x(t)} x2={x(t)} y1={10} y2={H - 18} stroke="var(--line-2)" /><text x={x(t)} y={H - 4} textAnchor="middle">{t.toFixed(2)}</text></g>)}
      <line x1={x(chance)} x2={x(chance)} y1={6} y2={H - 18} stroke="var(--ink)" strokeDasharray="2 3" />
      <text x={x(chance) + 4} y={12} style={{ fill: "var(--ink)" }}>chance {chance}</text>
      {dots.map((d) => <circle key={d.id} cx={d.cx} cy={y - d.level * 8} r={3.5} fill="var(--ink)" opacity={0.6}><title>sub-{d.id} {d.v.toFixed(3)}</title></circle>)}
      {Number.isFinite(st.lo) && <line x1={x(st.lo)} x2={x(st.hi)} y1={y + 10} y2={y + 10} stroke="var(--ink)" strokeWidth={2} />}
      <line x1={x(st.mean)} x2={x(st.mean)} y1={y + 4} y2={y + 16} stroke="var(--ink)" strokeWidth={2.5} />
      <text x={width - m.r + 12} y={y + 4} className="tnum">{st.mean.toFixed(3)} {Number.isFinite(st.lo) && <tspan style={{ fill: "var(--ink-3)" }}>95% CI {st.lo.toFixed(3)}–{st.hi.toFixed(3)}</tspan>}</text>
    </svg>
  );
}

export function GlobalView({ a, res, Provenance, used, excluded }) {
  const [all, setAll] = useState(false);
  const shown = all ? res.metrics : [res.metric];
  const chanceOf = (mt) => (mt === res.metric ? res.chance : mt === "AUC" ? 0 : null);
  return (
    <div className="page plain">
      <h1>{a.name}</h1>
      <p className="lede">{res.metric} over the whole mask, one value per subject. {res.chance_source ? <>Chance {res.chance}: {res.chance_source}. </> : null}<b>Descriptive only:</b> {(res.null?.description || "no test was run").replace(/^none:\s*/, "")}.</p>
      <Provenance a={a} res={res} used={used} excluded={excluded} />
      <section className="fig">
        <div className="fig-head"><h2>Every subject{all ? ", every metric" : ""}</h2>
          {res.metrics.length > 1 && <div className="seg"><button className={!all ? "on" : ""} onClick={() => setAll(false)}>{res.metric}</button><button className={all ? "on" : ""} onClick={() => setAll(true)}>All {res.metrics.length} metrics</button></div>}
        </div>
        {shown.map((mt) => {
          const rows = res.subjects.map((s) => ({ id: s.id, v: s.values[mt] })).filter((r) => typeof r.v === "number");
          return rows.length ? <MetricStrip key={mt} metric={mt} rows={rows} chance={chanceOf(mt) ?? Math.min(...rows.map((r) => r.v))} /> : null;
        })}
        <p className="fig-caption">One dot per subject (equal values stack); the tick is the mean and the bar its 95% t-interval over subjects. It describes the spread; it is not a test.</p>
      </section>
      <section className="fig">
        <div className="fig-head"><h2>The values</h2></div>
        <table className="values tnum">
          <thead><tr><th>subject</th>{res.metrics.map((mt) => <th key={mt}>{mt}</th>)}</tr></thead>
          <tbody>{res.subjects.map((s) => <tr key={s.id}><td>sub-{s.id}</td>{res.metrics.map((mt) => <td key={mt}>{typeof s.values[mt] === "number" ? s.values[mt].toFixed(3) : "n/a"}</td>)}</tr>)}</tbody>
        </table>
      </section>
    </div>
  );
}
