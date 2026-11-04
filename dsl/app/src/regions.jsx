// Regional results: parcel means on the cortex (surfview) or on slices of the
// atlas the analysis used, a network summary, and the selected parcel's
// subjects. Atlas semantics (region identity, hemisphere, projection of
// parcel values to voxels and vertices) come from neuroatlasjs.
import React, { Suspense, useEffect, useMemo, useState } from "react";
import { VolumeAtlas, SurfaceAtlas, ParcelValues } from "neuroatlasjs";
import { api, decodeU8 } from "./api.js";
import { ParcelSlices, DIVERGING } from "./brain.jsx";
import { meanCI, ticks } from "./charts.jsx";

const CortexViews = React.lazy(() => import("./surface.jsx"));

const b64 = (s) => { const bin = atob(s); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return u.buffer; };
const NETWORKS = [["Vis", "Visual"], ["SomMot", "Somatomotor"], ["DorsAttn", "Dorsal attention"], ["SalVentAttn", "Salience / ventral attention"],
  ["Limbic", "Limbic"], ["Cont", "Control"], ["Default", "Default"]];
const scale = (d0, d1, r0, r1) => (v) => r0 + ((v - d0) / (d1 - d0)) * (r1 - r0);

export function useRegionalData(res) {
  const [tpl, setTpl] = useState(null);
  const [surf, setSurf] = useState(null);
  useEffect(() => { api.template().then(setTpl); api.surface?.("fsaverage5").then(setSurf).catch(() => setSurf(false)); }, []);
  return useMemo(() => {
    if (!res?.regions || !tpl) return null;
    const dims = tpl.dims;
    const labels = new Int32Array(b64(res.labels));
    const regionDefs = res.regions.map((r) => ({ id: r.index, label: r.name, fullLabel: r.full, hemisphere: r.hemisphere, network: r.network }));
    // Studio template grid: 4 mm, origin from the fixture (see make-regions.R).
    const A = [[4, 0, 0, -96], [0, 4, 0, -132], [0, 0, 4, -78], [0, 0, 0, 1]];
    const space = { dim: dims, trans: { rows: 4, columns: 4, get: (r, c) => A[r][c] } };
    const volAtlas = new VolumeAtlas({ name: res.atlas.name, template: "MNI152", volume: { space, length: labels.length, getData: () => labels }, regions: regionDefs });
    const stats = res.regions.map((r) => ({ ...r, ...meanCI(r.values.filter((v) => typeof v === "number")) }));
    const pv = new ParcelValues(volAtlas.regions, stats.map((s) => ({ id: s.index, value: s.mean })));
    let surface = null;
    if (surf) {
      const hemi = (h) => ({
        positions: new Float32Array(b64(surf[h].positions)), faces: new Uint32Array(b64(surf[h].faces)),
        sulc: new Float32Array(b64(surf[h].sulc)), labels: new Int32Array(b64(surf[h].labels)),
      });
      const mesh = { left: hemi("left"), right: hemi("right") };
      const sAtlas = new SurfaceAtlas({ name: res.atlas.name, template: surf.template, regions: volAtlas.regions, left: mesh.left.labels, right: mesh.right.labels });
      const p = sAtlas.project(pv, NaN);
      const rel = (a) => Float32Array.from(a, (v) => v - res.chance);
      surface = { mesh, projected: { left: rel(p.left), right: rel(p.right) }, template: surf.template, source: surf.source };
    }
    return { dims, template: decodeU8(tpl.template), labels, voxelMeans: volAtlas.project(pv, NaN), stats, surface };
  }, [res, tpl, surf]);
}

// Each network is a row; each parcel a dot at its subject mean (LH filled,
// RH hollow). The selected parcel and its homologue are ringed.
function NetworkStrip({ stats, chance, selected, homologue, onSelect, hover, setHover, width = 1060 }) {
  const all = stats.map((s) => s.mean);
  const lo = Math.min(chance - 0.02, ...all) - 0.005, hi = Math.max(...all) + 0.01;
  const m = { l: 210, r: 24, t: 22, b: 26 }, rowH = 30;
  const H = m.t + NETWORKS.length * rowH + m.b;
  const x = scale(lo, hi, m.l, width - m.r);
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label="Parcel means by network">
      {ticks(lo, hi).map((t) => <g key={t}><line x1={x(t)} x2={x(t)} y1={m.t - 6} y2={H - m.b} stroke="var(--line-2)" /><text x={x(t)} y={H - 8} textAnchor="middle">{t.toFixed(2)}</text></g>)}
      <line x1={x(chance)} x2={x(chance)} y1={m.t - 10} y2={H - m.b} stroke="var(--ink)" strokeDasharray="2 3" />
      <text x={x(chance) + 4} y={m.t - 12} style={{ fill: "var(--ink)" }}>chance {chance}</text>
      {NETWORKS.map(([key, name], i) => {
        const rows = stats.filter((s) => s.network === key);
        return (
          <g key={key} transform={`translate(0, ${m.t + i * rowH})`}>
            <text x={m.l - 12} y={rowH / 2 + 4} textAnchor="end">{name} <tspan style={{ fill: "var(--ink-3)" }}>{rows.length}</tspan></text>
            {rows.map((s, j) => {
              const lh = s.hemisphere === "left";
              const on = s.index === selected || s.index === homologue || s.index === hover;
              return (
                <g key={s.index} onClick={() => onSelect(s.index)} onMouseEnter={() => setHover(s.index)} onMouseLeave={() => setHover(null)} style={{ cursor: "pointer" }}>
                  {on && <circle cx={x(s.mean)} cy={rowH / 2 + (lh ? -4 : 4)} r={7} fill="none" stroke="var(--ink)" strokeWidth={s.index === selected ? 2 : 1} />}
                  <circle cx={x(s.mean)} cy={rowH / 2 + (lh ? -4 : 4)} r={3} fill={lh ? "var(--ink)" : "var(--surface)"} stroke="var(--ink)" strokeWidth={lh ? 0 : 1} opacity={on ? 1 : 0.55} />
                </g>
              );
            })}
          </g>
        );
      })}
    </svg>
  );
}

// The selected parcel and its homologue: every subject, mean, 95% interval.
function ParcelDetail({ sel, hom, chance, subjects, width = 1060 }) {
  const rows = [sel, hom].filter(Boolean);
  const all = rows.flatMap((r) => r.values);
  const lo = Math.min(chance - 0.02, ...all) - 0.01, hi = Math.max(...all) + 0.01;
  const m = { l: 210, r: 230, t: 18, b: 26 }, rowH = 34;
  const H = m.t + rows.length * rowH + m.b;
  const x = scale(lo, hi, m.l, width - m.r);
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label="Selected parcel, every subject">
      {ticks(lo, hi).map((t) => <g key={t}><line x1={x(t)} x2={x(t)} y1={m.t - 6} y2={H - m.b} stroke="var(--line-2)" /><text x={x(t)} y={H - 8} textAnchor="middle">{t.toFixed(2)}</text></g>)}
      <line x1={x(chance)} x2={x(chance)} y1={m.t - 10} y2={H - m.b} stroke="var(--ink)" strokeDasharray="2 3" />
      {rows.map((r, i) => (
        <g key={r.index} transform={`translate(0, ${m.t + i * rowH})`}>
          <text x={m.l - 12} y={rowH / 2 + 4} textAnchor="end" style={{ fontWeight: i === 0 ? 600 : 400 }}>{r.hemisphere === "left" ? "LH" : "RH"} {r.name}</text>
          {r.values.map((v, j) => <circle key={j} cx={x(v)} cy={rowH / 2 + ((j * 7919) % 5 - 2) * 1.6} r={3} fill={r.hemisphere === "left" ? "var(--ink)" : "var(--surface)"} stroke="var(--ink)" strokeWidth={r.hemisphere === "left" ? 0 : 1} opacity={0.55}><title>sub-{subjects[j]} {v.toFixed(3)}</title></circle>)}
          <line x1={x(r.lo)} x2={x(r.hi)} y1={rowH / 2} y2={rowH / 2} stroke="var(--ink)" strokeWidth={2} />
          <line x1={x(r.mean)} x2={x(r.mean)} y1={rowH / 2 - 7} y2={rowH / 2 + 7} stroke="var(--ink)" strokeWidth={2.5} />
          <text x={width - m.r + 12} y={rowH / 2 + 4} className="tnum">{r.mean.toFixed(3)} <tspan style={{ fill: "var(--ink-3)" }}>95% CI {r.lo.toFixed(3)}–{r.hi.toFixed(3)}</tspan></text>
        </g>
      ))}
    </svg>
  );
}

function DivergingKey({ chance, span }) {
  const W = 260;
  const stops = DIVERGING.map((c, i) => `rgb(${c.join(",")}) ${(i / (DIVERGING.length - 1)) * 100}%`).join(",");
  const x = scale(chance - span, chance + span, 0, W);
  return (
    <div className="dkey">
      <div className="dkey-ramp" style={{ width: W, background: `linear-gradient(90deg, ${stops})` }} />
      <svg className="chart" width={W + 20} height={18}>
        {[chance - span, chance, chance + span].map((t) => <text key={t} x={Math.min(W - 2, Math.max(2, x(t)))} y={13} textAnchor={t === chance ? "middle" : t < chance ? "start" : "end"}>{t === chance ? `chance ${chance}` : t.toFixed(2)}</text>)}
      </svg>
    </div>
  );
}

export function RegionsView({ a, res, Provenance, used, excluded }) {
  const data = useRegionalData(res);
  const [view, setView] = useState("cortex");
  const [sel, setSel] = useState(null);
  const [hover, setHover] = useState(null);
  if (!data) return <div className="page">Loading…</div>;
  const { stats } = data;
  const byId = new Map(stats.map((s) => [s.index, s]));
  const top = [...stats].sort((p, q) => q.mean - p.mean)[0];
  const selId = sel ?? top.index;
  const selected = byId.get(selId);
  const homologue = stats.find((s) => s.name === selected.name && s.hemisphere !== selected.hemisphere);
  const span = 0.12;
  const pick = (id) => { if (id > 0 && byId.has(id)) setSel(id); };
  const pz = 18;
  const slices = [-6, -2, 2, 6, 10, 14].map((d) => pz + d);
  return (
    <div className="page plain">
      <h1>{a.name}</h1>
      <p className="lede">{res.metric} by parcel ({res.atlas.name}) · chance {res.chance} ({res.chance_source}) · <span className="stamp">descriptive only</span> rMVPA has no permutation path for transfer analyses.</p>
      <Provenance a={a} res={res} used={used} excluded={excluded} extra={<div><dt>atlas</dt><dd>Schaefer-200 · md5 {res.atlas.md5.slice(0, 8)}</dd></div>} />

      <section className="fig">
        <div className="fig-head">
          <h2>Mean transfer accuracy of {used} subjects, by parcel</h2>
          <div className="seg">
            <button className={view === "cortex" ? "on" : ""} onClick={() => setView("cortex")} disabled={!data.surface}>Cortex</button>
            <button className={view === "slices" ? "on" : ""} onClick={() => setView("slices")}>Slices</button>
          </div>
        </div>
        {view === "cortex" && data.surface ? (
          <Suspense fallback={<div className="cortex-wait">Drawing the cortex…</div>}>
            <CortexViews surface={data.surface.mesh} projected={data.surface.projected} selected={selId} span={span} onPick={pick} />
          </Suspense>
        ) : (
          <ParcelSlices template={data.template} values={data.voxelMeans} labels={data.labels} dims={data.dims} center={res.chance} span={span} slices={slices} selected={selId} onPick={pick} />
        )}
        <div className="subhead">
          <p className="fig-caption" style={{ margin: 0 }}>
            Each parcel takes one colour: the mean over subjects, blue below chance and red above. Click a parcel to see its subjects.
            {view === "cortex"
              ? <> Shown on {data.surface?.template} (pial, sulcal shading) with the atlas's surface labels; the analysis itself ran on the volumetric atlas in the data space. Medial wall uncoloured.</>
              : <> The volumetric atlas the analysis used, resampled to the studio template for display.</>}
          </p>
          <DivergingKey chance={res.chance} span={span} />
        </div>
      </section>

      <section className="fig">
        <div className="fig-head"><h2>{selected.hemisphere === "left" ? "LH" : "RH"} {selected.name} <span className="origin">{selected.network} network · parcel {selected.index}{sel == null ? " · opened on the highest mean" : ""}</span></h2></div>
        <ParcelDetail sel={selected} hom={homologue} chance={res.chance} subjects={res.subjects} />
        <p className="fig-caption">Every subject's accuracy in this parcel and its homologue; the tick is the mean and the bar its 95% t-interval. Opening on the highest mean is a selection: the interval for that parcel is optimistic.</p>
      </section>

      <section className="fig">
        <div className="fig-head"><h2>All {stats.length} parcels, by network</h2><span className="origin">LH ● · RH ○ · click a dot to select</span></div>
        <NetworkStrip stats={stats} chance={res.chance} selected={selId} homologue={homologue?.index} onSelect={setSel} hover={hover} setHover={setHover} />
        <p className="fig-caption">One dot per parcel at its subject mean. These are {stats.length} unadjusted means; some will sit above chance by luck.</p>
      </section>
    </div>
  );
}
