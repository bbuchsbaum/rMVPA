// Regional results: parcel means on the cortex (surfview) or on slices of the
// atlas the analysis used, a network summary, and the selected parcel's
// subjects. Atlas semantics (region identity, hemisphere, projection of
// parcel values to voxels and vertices) come from neuroatlasjs.
import React, { Suspense, useEffect, useMemo, useState } from "react";
import { VolumeAtlas, SurfaceAtlas, ParcelValues } from "neuroatlasjs";
import { api, decodeU8 } from "./api.js";
import { ParcelSlices, diverge } from "./brain.jsx";
import { meanCI, ticks } from "./charts.jsx";

const CortexViews = React.lazy(() => import("./surface.jsx"));

const b64 = (s) => { const bin = atob(s); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return u.buffer; };
const NETWORKS = [["Vis", "Visual"], ["SomMot", "Somatomotor"], ["DorsAttn", "Dorsal attention"], ["SalVentAttn", "Salience / ventral attention"],
  ["Limbic", "Limbic"], ["Cont", "Control"], ["Default", "Default"]];
const scale = (d0, d1, r0, r1) => (v) => r0 + ((v - d0) / (d1 - d0)) * (r1 - r0);

const FIXTURE_AFFINE = [[4, 0, 0, -96], [0, 4, 0, -132], [0, 0, 4, -78], [0, 0, 0, 1]]; // example template (see make-regions.R)
export const hemiTag = (s) => (s?.hemisphere === "left" ? "LH " : s?.hemisphere === "right" ? "RH " : "");

// Rows of the parcel strip: the atlas's networks when its label table names
// them (Schaefer's in their usual order first), else one row of all parcels.
function groupsOf(stats) {
  const nets = [...new Set(stats.map((s) => s.network).filter(Boolean))];
  if (!nets.length) return [[null, "All parcels"]];
  const known = NETWORKS.filter(([k]) => nets.includes(k));
  return [...known, ...nets.filter((k) => !known.some(([kk]) => kk === k)).map((k) => [k, k])];
}

export function useRegionalData(res) {
  const [tpl, setTpl] = useState(null);
  const [surf, setSurf] = useState(null);
  const surfName = res?.atlas?.surface;
  useEffect(() => { api.template().then(setTpl); }, []);
  useEffect(() => {
    setSurf(null);
    if (surfName) api.surface(surfName).then(setSurf).catch(() => setSurf(false)); else setSurf(false);
  }, [surfName]);
  return useMemo(() => {
    if (!res?.regions || !tpl) return null;
    const dims = tpl.dims;
    const labels = new Int32Array(b64(res.labels));
    // neuroatlasjs refuses empty descriptors, so leave out the ones the atlas lacks.
    const regionDefs = res.regions.map((r) => Object.fromEntries(Object.entries({ id: r.index, label: r.name, fullLabel: r.full, hemisphere: r.hemisphere, network: r.network }).filter(([, v]) => v != null && v !== "")));
    // Voxel-to-world affine of the template grid the labels were resampled to.
    const A = tpl.affine || FIXTURE_AFFINE;
    const space = { dim: dims, trans: { rows: 4, columns: 4, get: (r, c) => A[r][c] } };
    const volAtlas = new VolumeAtlas({ name: res.atlas.name, template: tpl.space || "MNI152", volume: { space, length: labels.length, getData: () => labels }, regions: regionDefs });
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
    // Six axial slices spread over the labelled extent.
    const plane = dims[0] * dims[1];
    let z0 = dims[2], z1 = -1;
    for (let i = 0; i < labels.length; i++) if (labels[i] > 0) { const z = Math.floor(i / plane); if (z < z0) z0 = z; if (z > z1) z1 = z; }
    if (z1 < 0) { z0 = 0; z1 = dims[2] - 1; }
    const slices = [...new Set(Array.from({ length: 6 }, (_, i) => Math.round(z0 + ((i + 1) / 7) * (z1 - z0))))];
    return { dims, template: decodeU8(tpl.template), labels, voxelMeans: volAtlas.project(pv, NaN), stats, surface, slices };
  }, [res, tpl, surf]);
}

// Each network is a row; each parcel a dot at its subject mean (LH filled,
// RH hollow). The selected parcel and its homologue are ringed.
function NetworkStrip({ stats, chance, selected, homologue, onSelect, hover, setHover, domain, width = 1060 }) {
  const [lo, hi] = domain;
  const m = { l: 210, r: 24, t: 22, b: 26 }, rowH = 30;
  const groups = groupsOf(stats);
  const H = m.t + groups.length * rowH + m.b;
  const x = scale(lo, hi, m.l, width - m.r);
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label="Parcel means by network">
      {ticks(lo, hi).map((t) => <g key={t}><line x1={x(t)} x2={x(t)} y1={m.t - 6} y2={H - m.b} stroke="var(--line-2)" /><text x={x(t)} y={H - 8} textAnchor="middle">{t.toFixed(2)}</text></g>)}
      <line x1={x(chance)} x2={x(chance)} y1={m.t - 10} y2={H - m.b} stroke="var(--ink)" strokeDasharray="2 3" />
      <text x={x(chance) + 4} y={m.t - 12} style={{ fill: "var(--ink)" }}>chance {chance}</text>
      {groups.map(([key, name], i) => {
        const rows = key == null ? stats : stats.filter((s) => s.network === key);
        const nm = rows.reduce((a, s) => a + s.mean, 0) / rows.length;
        const tagged = rows.filter((s) => s.index === selected || s.index === homologue);
        return (
          <g key={key} transform={`translate(0, ${m.t + i * rowH})`}>
            <text x={m.l - 12} y={rowH / 2 + 4} textAnchor="end">{name} <tspan style={{ fill: "var(--ink-3)" }}>{rows.length}</tspan></text>
            <line x1={x(nm)} x2={x(nm)} y1={3} y2={rowH - 3} stroke="var(--ink-3)" strokeWidth={2} />
            {tagged.map((s) => <text key={s.index} x={x(s.mean) + 10} y={rowH / 2 + (s.hemisphere === "left" ? -6 : 14)} style={{ fill: "var(--ink)", fontWeight: s.index === selected ? 600 : 400 }}>{hemiTag(s)}{s.name}</text>)}
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
function ParcelDetail({ sel, hom, chance, subjects, domain, width = 1060 }) {
  const rows = [sel, hom].filter(Boolean);
  const [lo, hi] = domain;
  const m = { l: 210, r: 230, t: 18, b: 26 }, rowH = 34;
  const H = m.t + rows.length * rowH + m.b;
  const x = scale(lo, hi, m.l, width - m.r);
  return (
    <svg className="chart" width={width} height={H} role="img" aria-label="Selected parcel, every subject">
      {ticks(lo, hi).map((t) => <g key={t}><line x1={x(t)} x2={x(t)} y1={m.t - 6} y2={H - m.b} stroke="var(--line-2)" /><text x={x(t)} y={H - 8} textAnchor="middle">{t.toFixed(2)}</text></g>)}
      <line x1={x(chance)} x2={x(chance)} y1={m.t - 10} y2={H - m.b} stroke="var(--ink)" strokeDasharray="2 3" />
      {rows.map((r, i) => (
        <g key={r.index} transform={`translate(0, ${m.t + i * rowH})`}>
          <text x={m.l - 12} y={rowH / 2 + 4} textAnchor="end" style={{ fontWeight: i === 0 ? 600 : 400 }}>{hemiTag(r)}{r.name}</text>
          {r.values.map((v, j) => <circle key={j} cx={x(v)} cy={rowH / 2} r={3} fill={r.hemisphere === "left" ? "var(--ink)" : "var(--surface)"} stroke="var(--ink)" strokeWidth={r.hemisphere === "left" ? 0 : 1} opacity={0.55}><title>sub-{subjects[j]} {v.toFixed(3)}</title></circle>)}
          <line x1={x(r.lo)} x2={x(r.hi)} y1={rowH / 2} y2={rowH / 2} stroke="var(--ink)" strokeWidth={2} />
          <line x1={x(r.mean)} x2={x(r.mean)} y1={rowH / 2 - 7} y2={rowH / 2 + 7} stroke="var(--ink)" strokeWidth={2.5} />
          <text x={width - m.r + 12} y={rowH / 2 + 4} className="tnum">{r.mean.toFixed(3)} <tspan style={{ fill: "var(--ink-3)" }}>95% CI {r.lo.toFixed(3)}–{r.hi.toFixed(3)}</tspan></text>
        </g>
      ))}
    </svg>
  );
}

function DivergingKey({ chance, span, deadzone }) {
  const W = 260;
  const stops = Array.from({ length: 41 }, (_, i) => `rgb(${diverge(i / 20 - 1, deadzone).join(",")}) ${i * 2.5}%`).join(",");
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
  const [viewPick, setView] = useState("cortex");
  const [sel, setSel] = useState(null);
  const [hover, setHover] = useState(null);
  if (!data) return <div className="page">Loading…</div>;
  const { stats } = data;
  const view = data.surface ? viewPick : "slices";
  const byId = new Map(stats.map((s) => [s.index, s]));
  const selected = sel != null ? byId.get(sel) : null;
  const homologue = selected?.hemisphere ? stats.find((s) => s.name === selected.name && s.hemisphere && s.hemisphere !== selected.hemisphere) : null;
  const marked = [selected?.index, homologue?.index].filter((v) => v != null);
  const names = Object.fromEntries(stats.map((s) => [s.index, `${hemiTag(s)}${s.name}`]));
  // Colour domain fitted to the parcel means; a ±0.01 band around chance stays neutral.
  const maxDev = Math.max(1e-6, ...stats.filter((s) => Number.isFinite(s.mean)).map((s) => Math.abs(s.mean - res.chance)));
  const span = Math.ceil(maxDev / 0.02) * 0.02;
  const deadzone = 0.01 / span;
  const pick = (id) => { if (id > 0 && byId.has(id)) setSel(id); };
  const slices = data.slices;
  const allVals = stats.flatMap((s) => s.values).filter((v) => typeof v === "number");
  const domain = [Math.min(res.chance, ...allVals) - 0.01, Math.max(...allVals) + 0.01];
  const above = stats.filter((s) => s.lo > res.chance).length;
  const ex = (p) => a.explain.find((r) => r.path === p)?.value;
  const direction = a.verb === "transfer" ? `, trained on ${ex("transfer.train_set")}, tested once on ${ex("transfer.test_set")}` : ` (${a.verb})`;
  const what = res.metric === "Accuracy" ? "accuracy" : res.metric;
  const tested = res.null && res.null.kind !== "none";
  return (
    <div className="page plain">
      <h1>{a.name}</h1>
      <p className="lede">{res.metric} by parcel{direction}. Atlas: {res.atlas.name}. {res.chance_source ? <>Chance {res.chance}: {res.chance_source}. </> : <>Reference {res.chance}. </>}{tested ? <>Per-parcel tests are in the region tables; this page describes the values.</> : <><b>Descriptive only:</b> {(res.null?.description || "no test was run").replace(/^none:\s*/, "")}.</>}</p>
      <Provenance a={a} res={res} used={used} excluded={excluded} extra={<div><dt>atlas</dt><dd>{res.atlas.short || res.atlas.name}{res.atlas.md5 ? ` · md5 ${res.atlas.md5.slice(0, 8)}` : ""}</dd></div>} />

      <section className="fig">
        <div className="fig-head">
          <h2>Mean {what} of {res.subjects.length} subjects in each parcel</h2>
          <div className="seg">
            <button className={view === "cortex" ? "on" : ""} onClick={() => setView("cortex")} disabled={!data.surface}>Cortex</button>
            <button className={view === "slices" ? "on" : ""} onClick={() => setView("slices")}>Slices</button>
          </div>
        </div>
        {view === "cortex" && data.surface ? (
          <Suspense fallback={<div className="cortex-wait">Drawing the cortex…</div>}>
            <CortexViews surface={data.surface.mesh} projected={data.surface.projected} selected={marked} span={span} deadzone={deadzone} onPick={pick} />
          </Suspense>
        ) : (
          <ParcelSlices template={data.template} values={data.voxelMeans} labels={data.labels} dims={data.dims} center={res.chance} span={span} deadzone={deadzone} slices={slices} selected={marked} names={names} onPick={pick} />
        )}
        <div className="subhead">
          <p className="fig-caption" style={{ margin: 0 }}>
            Each parcel takes one colour: the mean over subjects, blue below chance and red above; within ±0.01 of chance stays neutral. Click a parcel to see its subjects and its homologue (both outlined).
            {view === "cortex"
              ? <> Shown on {data.surface?.template} (pial, sulcal shading) with the atlas's surface labels; the analysis itself ran on the volumetric atlas in the data space. Medial wall uncoloured.</>
              : <> The volumetric atlas the analysis used, resampled to the studio template for display.</>}
          </p>
          <DivergingKey chance={res.chance} span={span} deadzone={deadzone} />
        </div>
      </section>

      <section className="fig">
        {selected ? <>
          <div className="fig-head"><h2>{names[selected.index]} <span className="origin">{selected.network ? `${selected.network} network · ` : ""}parcel {selected.index}</span></h2><button className="fix" onClick={() => setSel(null)}>Clear</button></div>
          <ParcelDetail sel={selected} hom={homologue} chance={res.chance} subjects={res.subjects} domain={domain} />
          <p className="fig-caption">Every subject's {what} in this parcel and its homologue; the tick is the mean and the bar its 95% t-interval. If you chose this parcel because it looked strong, its interval is optimistic.</p>
        </> : <div className="empty-detail"><b>No parcel selected.</b> Click a parcel on the brain or a dot below to see every subject's {what} there and in its homologue. The studio does not open on the strongest parcel, so the first one you see is your choice.</div>}
      </section>

      <section className="fig">
        <div className="fig-head"><h2>All {stats.length} parcels{stats.some((s) => s.network) ? ", by network" : ""}</h2><span className="origin">LH ● · RH ○ · grey tick: {stats.some((s) => s.network) ? "network" : "overall"} mean · click a dot to select</span></div>
        <NetworkStrip stats={stats} chance={res.chance} selected={selected?.index} homologue={homologue?.index} onSelect={setSel} hover={hover} setHover={setHover} domain={domain} />
        <p className="fig-caption"><b>{above} of {stats.length}</b> parcels have a 95% interval entirely above chance. With per-parcel intervals and no effect anywhere, about {Math.round(stats.length * 0.025)} would by chance. This count is a description, not a test. Same axis as the panel above.</p>
      </section>
    </div>
  );
}
