import React, { useEffect, useMemo, useState } from "react";
import { api, isExample, decodeF32, decodeU8 } from "./api.js";
import { Slices, PlaneSlice } from "./brain.jsx";
import { SharedAxis, RegionDots } from "./charts.jsx";

const VERBS = {
  decode: "Predict a label or value", transfer: "Train on one set, test on another",
  rsa: "Compare patterns to model RDMs", rsa_vector: "One RDM, trial vectors",
  rsa_contrast: "Contrast RSA (MS-ReVE)", rsa_feature: "Predict a feature space",
  rsa_era: "Encoding–retrieval similarity", rsa_connect: "Representational connectivity",
  rsa_mediate: "Representational mediation", encode: "Banded-ridge encoding",
  pattern: "Reduced-rank pattern model", manova: "Multivariate test of terms",
};

function useHash() {
  const [h, setH] = useState(location.hash.slice(1));
  useEffect(() => {
    const f = () => setH(location.hash.slice(1));
    addEventListener("hashchange", f);
    return () => removeEventListener("hashchange", f);
  }, []);
  return h;
}
const go = (h) => { location.hash = h; };
const groupBy = (xs, f) => xs.reduce((o, x) => ((o[f(x)] ||= []).push(x), o), {});

// The check problems of an analysis, grouped by message; and the YAML path a
// message points at, so the spec panel can mark the line.
function problems(project, name) {
  const errs = project.check.filter((c) => c.analysis === name && c.level !== "ok");
  return Object.values(groupBy(errs, (e) => e.message)).map((g) => {
    const m = g[0].message.match(/^assert (\w+):/);
    return { level: g[0].level, message: g[0].message, subjects: g.length, path: m ? "folds.assert" : null, assertion: m?.[1] };
  });
}

export default function App() {
  const [project, setProject] = useState(null);
  const hash = useHash();
  useEffect(() => { api.project().then(setProject); }, []);
  if (!project) return <div className="app"><div className="header"><div className="wordmark">rmvpa <span>studio</span></div></div></div>;

  const [kind, name] = hash.includes("-") ? [hash.slice(0, hash.indexOf("-")), hash.slice(hash.indexOf("-") + 1)] : ["project", null];
  const analysis = project.analyses.find((a) => a.name === name);
  const blocked = project.analyses.filter((a) => problems(project, a.name).some((p) => p.level === "error"));

  return (
    <div className="app">
      <header className="header">
        <div className="wordmark">rmvpa <span>studio</span></div>
        <div className="crumbs">
          <b>{project.name}</b>
          <span className="mono">{project.root}</span>
          {isExample && <span className="example-tag">example data</span>}
        </div>
        <button className={`pill ${blocked.length ? "bad" : ""}`} onClick={() => blocked.length && go(`a-${blocked[0].name}`)}><i />
          {blocked.length ? `${blocked.map((b) => b.name).join(", ")} blocked by check` : "Check passed"}
        </button>
        <button className="button primary">Run<kbd style={{ borderColor: "#555", color: "#ccc" }}>⌘↵</kbd></button>
      </header>
      <div className="body">
        <Nav project={project} current={name} kind={kind} blocked={blocked.map((b) => b.name)} />
        <main className="main">
          {kind === "a" && analysis ? <AnalysisView project={project} a={analysis} />
            : kind === "r" && analysis ? <ResultsView project={project} a={analysis} />
            : <ProjectView project={project} />}
        </main>
      </div>
    </div>
  );
}

function Nav({ project, current, kind, blocked }) {
  return (
    <nav className="nav">
      <h3>Project</h3>
      <button className={`nav-item ${kind === "project" ? "on" : ""}`} onClick={() => go("")}>Overview</button>
      <h3>Data sets</h3>
      {project.sets.map((s) => (
        <div key={s.name} className="nav-item">{s.name}<span className="sub">{s.subjects} subj · {s.runs} runs</span></div>
      ))}
      <h3>Analyses</h3>
      {project.analyses.map((a) => (
        <button key={a.name} className={`nav-item ${current === a.name ? "on" : ""} ${blocked.includes(a.name) ? "blocked" : ""}`} onClick={() => go(`a-${a.name}`)}>
          <span>{a.name}</span>
          <span className="dots" aria-hidden>{a.cells.map((c, i) => <i key={i} className={blocked.includes(a.name) ? "blocked" : c.status} />)}</span>
        </button>
      ))}
      <button className="nav-item" style={{ color: "var(--ink-3)" }}>+ New analysis<kbd>N</kbd></button>
    </nav>
  );
}

function Legend() {
  return (
    <div className="legend">
      <span><i className="l-done" />done</span>
      <span><i className="l-stale" />stale: spec or inputs changed</span>
      <span><i className="l-missing" />not run</span>
      <span><i className="l-blocked" />blocked by check</span>
    </div>
  );
}

function ProjectView({ project }) {
  const all = project.analyses.flatMap((a) => problems(project, a.name).map((p) => ({ ...p, analysis: a.name })));
  return (
    <div className="page">
      <h1>{project.name}</h1>
      <p className="lede">{project.analyses.length} analyses · {project.subjects.length} subjects · {project.sets.length} data sets. Each cell is one subject in one analysis.</p>
      <div className="stack">
        <div className="card">
          <table className="matrix">
            <thead><tr><th />{project.subjects.map((s) => <th key={s} className="tnum">{s}</th>)}</tr></thead>
            <tbody>
              {project.analyses.map((a) => {
                const isBlocked = all.some((p) => p.analysis === a.name && p.level === "error");
                return (
                  <tr key={a.name}>
                    <td className="name"><a href={`#a-${a.name}`}>{a.name}</a> <span className="verb">{a.verb}</span></td>
                    {a.cells.map((c) => <td key={c.subject}><span className={`cell ${isBlocked ? "blocked" : c.status}`} title={`sub-${c.subject}: ${isBlocked ? "blocked" : c.status}`} /></td>)}
                  </tr>
                );
              })}
            </tbody>
          </table>
          <Legend />
        </div>
        <div className="card">
          <h2>Check</h2>
          {all.length === 0 ? <p>No problems.</p> : all.map((p) => (
            <div className="finding" key={p.analysis + p.message}>
              <span className={`lvl ${p.level}`}>{p.level.toUpperCase()}</span>
              <a className="mono" href={`#a-${p.analysis}`}>{p.analysis}</a>
              <span>{p.message} <span style={{ color: "var(--ink-3)" }}>· all {p.subjects} subjects</span></span>
              {p.assertion && <button className="fix">Remove assert</button>}
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

function AnalysisView({ project, a }) {
  const get = (p) => a.explain.find((r) => r.path === p);
  const probs = problems(project, a.name);
  const blocking = probs.find((p) => p.level === "error");
  const inf = a.inference;
  const hasResults = !blocking && a.cells.some((c) => c.status === "done");
  const verbFields = a.explain.filter((r) => r.depth === 1 && r.path.startsWith(a.verb + "."));
  const Field = ({ r }) => r ? (
    <div className="field">
      <label htmlFor={r.path}>{r.key}</label>
      <select id={r.path} className="select" defaultValue={r.value}><option>{r.value}</option></select>
      <span className={`origin ${(r.origin || "").split("@")[0].split(" ")[0]}`}>{r.origin}</span>
    </div>
  ) : null;
  return (
    <div className="analysis">
      <div className="spine">
        <h1 style={{ margin: "0 0 4px", fontSize: 22 }}>{a.name}</h1>
        <p className="lede">analyses/{a.name}.yml · {hasResults ? <a href={`#r-${a.name}`}>Results →</a> : <a className="disabled">No results yet</a>}</p>
        {blocking && (
          <div className="banner" role="alert">
            <div className="grow"><b>Blocked by check.</b> <code>assert {blocking.assertion}</code> does not apply to {a.verb}. Nothing will run until it is fixed.</div>
            <button className="fix">Remove assert</button>
            <button className="fix">Open in YAML</button>
          </div>
        )}
        <Section n={1} title="Question">
          <div className="qline"><b>{a.verb}</b><span>{VERBS[a.verb]}</span><button className="fix">Change<kbd>⌘K</kbd></button></div>
        </Section>
        <Section n={2} title="Data">
          {verbFields.map((r) => <Field key={r.path} r={r} />)}
        </Section>
        <Section n={3} title="Folds" badge={blocking?.path === "folds.assert" ? "1 problem" : null}>
          <Field r={get("folds.by")} /><Field r={get("folds.scheme")} /><Field r={get("folds.assert")} />
        </Section>
        <Section n={4} title="Scope">
          <div className="field"><label>kind</label><select className="select"><option>{a.scope}</option></select><span className="origin analysis">analysis</span></div>
          {a.explain.filter((r) => r.path.startsWith(`scope.${a.scope}.`)).map((r) => <Field key={r.path} r={r} />)}
        </Section>
        <Section n={5} title="Inference">
          <p className={`note ${inf.kind === "none" ? "warn" : ""}`}>{inf.description}</p>
          {inf.kind !== "none" && <><Field r={get("inference.permutations")} /><Field r={get("inference.correction")} /></>}
          {inf.kind === "searchlight_permutation" && inf.permutations < 1000 && (
            <p className="advice">{inf.permutations} permutations pool with about 10% of centres into a null of ~250,000 values, but neighbouring centres are correlated, so the null behaves more like {inf.permutations} independent draws. Use at least 1,000 before reporting FDR-corrected maps.</p>
          )}
        </Section>
      </div>
      <aside className="spec">
        <h3>analyses/{a.name}.yml · resolved</h3>
        <div className="yaml">
          <div><span><span className="k">extends: </span><span className="v">../rmvpa.yml</span></span><span /></div>
          {a.explain.map((r) => {
            const dflt = r.origin?.startsWith("default");
            const err = blocking && r.path === blocking.path;
            return (
              <div key={r.path} className={err ? "err" : ""} title={r.value || ""}>
                <span>
                  <span className={dflt ? "dflt" : "k"}>{"  ".repeat(r.depth)}{r.key}:</span>
                  {r.value != null && <span className={dflt ? "dflt" : "v"}> {r.value}</span>}
                </span>
                <span className={`o ${(r.origin || "").split("@")[0].split(" ")[0]}`}>{r.origin}</span>
              </div>
            );
          })}
        </div>
      </aside>
    </div>
  );
}

function Section({ n, title, badge, children }) {
  return (
    <section className="section">
      <div className="rail"><span className="num">{n}</span></div>
      <div className="sec-body">
        <div className="sec-head"><h2>{title}</h2>{badge && <span className="badge">{badge}</span>}</div>
        {children}
      </div>
    </section>
  );
}

function Provenance({ a, res, used, excluded, include, setInclude }) {
  const sc = res.sidecar || {};
  return (
    <dl className="prov">
      <div><dt>subjects</dt><dd>n = {used}</dd></div>
      {excluded.length > 0 && <div><dt>excluded</dt><dd>{excluded.join(", ")}</dd></div>}
      {setInclude && <label className="toggle"><input id="include-stale" type="checkbox" checked={include} onChange={(e) => setInclude(e.target.checked)} /> include stale</label>}
      <div><dt>lock</dt><dd className={`lock-${a.lock}`}>{a.lock}</dd></div>
      {sc.spec_hash && <div><dt>spec</dt><dd>{sc.spec_hash}</dd></div>}
      {sc.inference && <div><dt>null</dt><dd>{sc.inference.permutations} perms, {sc.inference.shuffle}</dd></div>}
      <div><dt>rmvpa</dt><dd>{sc.rmvpa_version || "0.1.3"}</dd></div>
      <div><a href="#">All sidecars</a></div>
    </dl>
  );
}

function ResultsView({ project, a }) {
  const [res, setRes] = useState(null);
  const [tpl, setTpl] = useState(null);
  const [hover, setHover] = useState(null);
  const [include, setInclude] = useState(false);
  const [order, setOrder] = useState("atlas");
  useEffect(() => { setRes(null); api.results(a.name).then(setRes).catch(() => setRes({ missing: true })); api.template().then(setTpl); }, [a.name]);
  const vol = useMemo(() => (tpl ? { dims: tpl.dims, template: decodeU8(tpl.template) } : null), [tpl]);
  const map = useMemo(() => (res?.maps ? decodeF32(Object.values(res.maps)[0]) : null), [res]);
  const planes = useMemo(() => Object.fromEntries((res?.subjects || []).filter((s) => typeof s.slice === "string").map((s) => [s.id, decodeF32(s.slice)])), [res]);
  if (!res || !vol) return <div className="page">Loading…</div>;
  if (res.missing) return <div className="page"><h1>{a.name}</h1><p className="lede">No results yet.</p></div>;

  const statusOf = Object.fromEntries(a.cells.map((c) => [c.subject, c.status]));
  const excluded = a.cells.filter((c) => c.status === "missing" || (c.status === "stale" && !include)).map((c) => `sub-${c.subject} (${c.status === "missing" ? "not run" : "stale"})`);
  const used = a.cells.filter((c) => c.status === "done" || (c.status === "stale" && include)).length;

  if (res.regions) {
    return (
      <div className="page">
        <h1>{a.name}</h1>
        <p className="lede">{res.metric} by region · chance {res.chance}</p>
        <Provenance a={a} res={res} used={used} excluded={excluded} />
        <div className="card">
          <div className="fig-head">
            <h2>Transfer accuracy, every subject <span className="stamp" style={{ marginLeft: 8 }}>descriptive only</span></h2>
            <div className="seg">{["atlas", "mean"].map((o) => <button key={o} className={order === o ? "on" : ""} onClick={() => setOrder(o)}>{o === "atlas" ? "Atlas order" : "By mean"}</button>)}</div>
          </div>
          <RegionDots regions={res.regions} chance={res.chance} order={order} />
          <p className="fig-caption">rMVPA has no permutation path for transfer analyses, so these are descriptive: a dot per subject, the tick is the subject mean. Regions are listed in atlas order so the order carries no selection.</p>
        </div>
      </div>
    );
  }

  const fdr = res.fdr.threshold;
  const pz = res.peak_z;
  const slices = [-3, 0, 3, 6, 9, 12].map((d) => pz + d);
  const subjects = res.subjects.map((s) => ({ ...s, status: statusOf[s.id] || s.status }));
  return (
    <div className="page">
      <h1>{a.name}</h1>
      <p className="lede">{res.metric} · {res.classes} classes, chance {res.chance} · {res.null.description}</p>
      <Provenance a={a} res={res} used={used} excluded={excluded} include={include} setInclude={setInclude} />
      <div className="stack">
        <div className="card">
          <div className="fig-head"><h2>Group mean accuracy, slices around the peak</h2></div>
          <Slices template={vol.template} map={map} dims={vol.dims} threshold={fdr} top={fdr + 0.07} floor={res.chance} slices={slices} />
          <p className="fig-caption">
            <span className="swatch" style={{ background: "#d4492c" }} />FDR q &lt; .05 against the pooled null: accuracy ≥ {fdr.toFixed(3)}, {res.fdr.surviving} of {res.centres.n.toLocaleString()} centres.{" "}
            <span className="swatch" style={{ background: "#e8b9ac" }} />Faint: above chance, not significant. Smallest attainable p {res.null.p_floor}.
          </p>
        </div>
        <div className="card">
          <div className="fig-head"><h2>Where the map's threshold comes from, and every subject</h2></div>
          <SharedAxis res={{ ...res, subjects }} include={include} hover={hover} setHover={setHover} />
          <p className="fig-caption">Each subject's value is its accuracy in the region found from the <i>other</i> subjects' group map (leave-one-subject-out), so no subject chooses its own voxels. {res.group.note}.</p>
        </div>
        <div className="card">
          <div className="fig-head"><h2>Each subject at z {pz}</h2><span className="origin">same threshold as the group map</span></div>
          <div className="multiples">
            {a.cells.map((c) => {
              const st = statusOf[c.subject];
              if (!planes[c.subject]) return <div key={c.subject} className="tile missing">sub-{c.subject}<br />not run</div>;
              return (
                <div key={c.subject} className={`tile ${st === "stale" ? "stale" : ""} ${hover === c.subject ? "hl" : ""}`}
                  onMouseEnter={() => setHover(c.subject)} onMouseLeave={() => setHover(null)}>
                  <PlaneSlice template={vol.template} plane={planes[c.subject]} dims={vol.dims} k={pz} threshold={fdr} top={fdr + 0.07} floor={res.chance} />
                  <div className="lab"><span>{c.subject}</span><span>{st === "stale" ? "stale" : ""}</span></div>
                </div>
              );
            })}
          </div>
        </div>
      </div>
    </div>
  );
}
