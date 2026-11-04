import React, { useEffect, useMemo, useState } from "react";
import { api, isExample, decodeF32, decodeU8 } from "./api.js";
import { Slices, PlaneSlice, COUNT_BINS } from "./brain.jsx";
import { SubjectNull, SubjectStrip, AccuracyKey, RegionPairs, CountKey } from "./charts.jsx";

const VERBS = {
  decode: "Predict a label or value", transfer: "Train on one set, test on another",
  rsa: "Compare patterns to model RDMs", rsa_vector: "One RDM, trial vectors",
  rsa_contrast: "Contrast RSA (MS-ReVE)", rsa_feature: "Predict a feature space",
  rsa_era: "Encoding–retrieval similarity", rsa_connect: "Representational connectivity",
  rsa_mediate: "Representational mediation", encode: "Banded-ridge encoding",
  pattern: "Reduced-rank pattern model", manova: "Multivariate test of terms",
};
// The fold assertion that does apply, per verb, when a project-wide one doesn't.
const GUARD = { rsa_era: "train_test_disjoint_runs", transfer: "train_test_disjoint_runs" };
const NUMERIC = new Set(["radius", "iterations", "permutations", "subsample", "workers", "seed", "nfolds", "repeats", "max_rank", "rank"]);
const UNITS = { radius: "mm" };

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

// Check findings of an analysis, grouped by message, with the YAML path and
// the fix each one offers.
function problems(project, name) {
  const a = project.analyses.find((x) => x.name === name);
  const errs = project.check.filter((c) => c.analysis === name && c.level !== "ok");
  return Object.values(groupBy(errs, (e) => e.message)).map((g) => {
    const msg = g[0].message;
    const asrt = msg.match(/^assert (\w+):/);
    const perm = msg.match(/^inference\.permutations = (\d+)/);
    return {
      level: g[0].level, message: msg, subjects: g.length,
      path: asrt ? "folds.assert" : perm ? "inference.permutations" : null,
      fix: asrt ? (GUARD[a.verb] ? `Use ${GUARD[a.verb]} instead` : "Remove assert") : perm ? "Set permutations: 1000" : null,
      assertion: asrt?.[1], guard: GUARD[a.verb],
    };
  });
}

export default function App() {
  const [project, setProject] = useState(null);
  const [runHover, setRunHover] = useState(false);
  const hash = useHash();
  useEffect(() => { api.project().then(setProject); }, []);
  if (!project) return <div className="app"><div className="header"><div className="wordmark">rmvpa <span>studio</span></div></div></div>;

  const [kind, name] = hash.includes("-") ? [hash.slice(0, hash.indexOf("-")), hash.slice(hash.indexOf("-") + 1)] : ["project", null];
  const analysis = project.analyses.find((a) => a.name === name);
  const blocked = project.analyses.filter((a) => problems(project, a.name).some((p) => p.level === "error")).map((a) => a.name);
  const todoBy = project.analyses.filter((a) => !blocked.includes(a.name))
    .map((a) => [a.name, a.cells.filter((c) => c.status !== "done").length]).filter(([, n]) => n > 0);
  const todo = todoBy.reduce((s, [, n]) => s + n, 0);

  return (
    <div className="app">
      <header className="header">
        <div className="wordmark">rmvpa <span>studio</span></div>
        <div className="crumbs">
          <b>{project.name}</b>
          <span className="mono">{project.root}</span>
          {isExample && <span className="example-tag">example data</span>}
        </div>
        <button className={`pill ${blocked.length ? "bad" : ""}`} onClick={() => blocked.length && go(`a-${blocked[0]}`)}><i />
          {blocked.length ? `${blocked.join(", ")} blocked by check` : "Check passed"}
        </button>
        <button className="button primary" title={todoBy.map(([n, k]) => `${k} × ${n}`).join("\n")}
          onMouseEnter={() => { setRunHover(true); if (kind !== "project") go(""); }} onMouseLeave={() => setRunHover(false)} onFocus={() => setRunHover(true)} onBlur={() => setRunHover(false)}>
          Run {todo} cells{blocked.length ? ` · skips ${blocked.join(", ")}` : ""}<kbd style={{ borderColor: "#555", color: "#ccc" }}>⌘↵</kbd>
        </button>
      </header>
      <div className="body">
        <Nav project={project} current={name} kind={kind} blocked={blocked} />
        <main className="main">
          {kind === "a" && analysis ? <AnalysisView project={project} a={analysis} />
            : kind === "r" && analysis ? <ResultsView project={project} a={analysis} />
            : <ProjectView project={project} blocked={blocked} runHover={runHover} />}
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

function ProjectView({ project, blocked, runHover }) {
  const all = project.analyses.flatMap((a) => problems(project, a.name).map((p) => ({ ...p, analysis: a.name })));
  return (
    <div className="page plain">
      <h1>{project.name}</h1>
      <p className="lede">{project.analyses.length} analyses · {project.subjects.length} subjects · {project.sets.length} data sets. Each cell is one subject in one analysis.</p>
      <div className="stack">
        <div className="card">
          <table className="matrix">
            <thead><tr><th />{project.subjects.map((s) => <th key={s} className="tnum">{s}</th>)}<th /></tr></thead>
            <tbody>
              {project.analyses.map((a) => {
                const isBlocked = blocked.includes(a.name);
                const done = a.cells.filter((c) => c.status === "done").length;
                return (
                  <tr key={a.name}>
                    <td className="name"><a href={`#a-${a.name}`}>{a.name}</a> <span className="verb">{a.verb}</span></td>
                    {a.cells.map((c) => <td key={c.subject}><span className={`cell ${isBlocked ? "blocked" : c.status} ${runHover && !isBlocked && c.status !== "done" ? "target" : ""}`} title={`sub-${c.subject}: ${isBlocked ? "blocked" : c.status}`} /></td>)}
                    <td className="count tnum">{isBlocked ? "blocked" : `${done}/${a.cells.length}`}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
          <Legend />
          {runHover && <p className="fig-caption">Outlined: the cells Run will compute. Blocked analyses are skipped.</p>}
        </div>
        <div className="card">
          <h2>Check</h2>
          {all.length === 0 ? <p>No problems.</p> : all.map((p) => (
            <div className="finding" key={p.analysis + p.message}>
              <span className={`lvl ${p.level}`}>{p.level === "error" ? "Blocks" : "Advice"}</span>
              <a className="mono" href={`#a-${p.analysis}`}>{p.analysis}</a>
              <span>{p.message}</span>
              {p.fix && <button className="fix">{p.fix}</button>}
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
  const advice = probs.filter((p) => p.level === "warning");
  const inf = a.inference;
  const hasResults = !blocking && a.cells.some((c) => c.status === "done");
  const verbFields = a.explain.filter((r) => r.depth === 1 && r.path.startsWith(a.verb + "."));
  const errPaths = new Set(probs.map((p) => p.path).filter(Boolean));
  const Field = ({ r }) => r ? (
    <div className={`field ${blocking && blocking.path === r.path ? "err" : ""}`}>
      <label htmlFor={r.path}>{r.key}</label>
      {NUMERIC.has(r.key)
        ? <span className="stepper"><button type="button" aria-label={`decrease ${r.key}`}>−</button><input id={r.path} className="tnum" type="number" defaultValue={r.value} /><button type="button" aria-label={`increase ${r.key}`}>+</button>{UNITS[r.key] && <span className="unit">{UNITS[r.key]}</span>}</span>
        : <select id={r.path} className="select" defaultValue={r.value}><option>{r.value}</option></select>}
      <span className={`origin ${(r.origin || "").split("@")[0].split(" ")[0]}`}>{r.origin}</span>
    </div>
  ) : null;
  return (
    <div className="analysis-layout">
      <div className="spine">
        <h1 style={{ margin: "0 0 4px", fontSize: 22 }}>{a.name}</h1>
        <p className="lede">analyses/{a.name}.yml · {hasResults ? <a href={`#r-${a.name}`}>Results →</a> : <a className="disabled">No results yet</a>}</p>
        {blocking && (
          <div className="banner" role="alert">
            <div><b>Blocked by check.</b> <code>assert {blocking.assertion}</code> does not apply to {a.verb}.
              {blocking.guard ? <> The guard that does is <code>{blocking.guard}</code>: train and test sets must not share runs.</> : " No assertion of this kind applies."}
              {" "}The change is written to the analysis file and stamped in each sidecar.</div>
            <div className="banner-actions"><button className="fix strong">{blocking.fix}</button><button className="fix">Open in YAML</button></div>
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
        <Section n={5} title="Inference" badge={advice.length ? "advice" : null}>
          <p className="note">{inf.description}</p>
          {inf.kind !== "none" && <><Field r={get("inference.permutations")} /><Field r={get("inference.correction")} /></>}
          {advice.map((p) => <p key={p.message} className="advice">{p.message}. <button className="fix" style={{ marginLeft: 6 }}>{p.fix}</button></p>)}
        </Section>
      </div>
      <aside className="spec">
        <h3>analyses/{a.name}.yml · resolved</h3>
        <div className="yaml">
          <div><span><span className="k">extends: </span><span className="v">../rmvpa.yml</span></span><span /></div>
          {a.explain.map((r) => {
            const dflt = r.origin?.startsWith("default");
            const err = errPaths.has(r.path);
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

function Provenance({ a, res, used, excluded, extra }) {
  const sc = res.sidecar || {};
  return (
    <dl className="prov">
      <div><dt>n</dt><dd>{used}</dd></div>
      {excluded.length > 0 && <div><dt>excluded</dt><dd>{excluded.join(", ")}</dd></div>}
      <div><dt>lock</dt><dd className={`lock-${a.lock}`}>{a.lock === "none" ? "none · not preregistered" : a.lock}</dd></div>
      {sc.spec_hash && <div><dt>spec</dt><dd>{sc.spec_hash.slice(0, 10)}</dd></div>}
      {sc.inference && <div><dt>null</dt><dd>{sc.inference.permutations} perms · {sc.inference.shuffle}</dd></div>}
      {extra}
      <div><dt>rMVPA</dt><dd>{sc.rmvpa_version || "0.1.3"}</dd></div>
      <div><a href="#">Sidecars ↗</a></div>
    </dl>
  );
}

function ResultsView({ project, a }) {
  const [res, setRes] = useState(null);
  const [tpl, setTpl] = useState(null);
  const [hover, setHover] = useState(null);
  const [include, setInclude] = useState(false);
  const [order, setOrder] = useState("atlas");
  const [view, setView] = useState("count");
  const [tileOrder, setTileOrder] = useState("id");
  const [sel, setSel] = useState(null);
  useEffect(() => { setRes(null); api.results(a.name).then(setRes).catch(() => setRes({ missing: true })); api.template().then(setTpl); }, [a.name]);
  const vol = useMemo(() => (tpl ? { dims: tpl.dims, template: decodeU8(tpl.template) } : null), [tpl]);
  const maps = useMemo(() => (res?.maps ? Object.fromEntries(Object.entries(res.maps).map(([k, v]) => [k, decodeF32(v)])) : {}), [res]);
  const planes = useMemo(() => Object.fromEntries((res?.subjects || []).filter((s) => typeof s.slice === "string").map((s) => [s.id, decodeF32(s.slice)])), [res]);
  const statusOf = Object.fromEntries(a.cells.map((c) => [c.subject, c.status]));
  const ran = (res?.subjects || []).filter((s) => typeof s.slice === "string").map((s) => s.id);
  // ←/→ step through subjects.
  useEffect(() => {
    const f = (e) => {
      if (!ran.length || !["ArrowLeft", "ArrowRight"].includes(e.key) || e.target.closest("input,select")) return;
      const i = ran.indexOf(sel);
      setSel(ran[(i + (e.key === "ArrowRight" ? 1 : ran.length - 1)) % ran.length]);
    };
    addEventListener("keydown", f);
    return () => removeEventListener("keydown", f);
  }, [sel, ran.join()]);
  if (!res || !vol) return <div className="page">Loading…</div>;
  if (res.missing) return <div className="page"><h1>{a.name}</h1><p className="lede">No results yet.</p></div>;

  // One n for every panel: done cells, plus stale ones when included.
  const inN = (id) => statusOf[id] === "done" || (statusOf[id] === "stale" && include);
  const excluded = a.cells.filter((c) => !inN(c.subject)).map((c) => `sub-${c.subject} (${c.status === "missing" ? "not run" : c.status})`);
  const n = a.cells.filter((c) => inN(c.subject)).length;

  if (res.regions) {
    const intervals = res.regions.length;
    return (
      <div className="page plain">
        <h1>{a.name}</h1>
        <p className="lede">{res.metric} by region · chance {res.chance} ({res.chance_source}) · <span className="stamp">descriptive only</span> rMVPA has no permutation path for transfer analyses.</p>
        <Provenance a={a} res={res} used={n} excluded={excluded} extra={<div><dt>atlas</dt><dd>{res.atlas.file.split("/").pop()} · md5 {res.atlas.md5.slice(0, 8)}</dd></div>} />
        <section className="fig">
          <div className="fig-head">
            <h2>Transfer accuracy, every subject, hemispheres paired
              {order === "mean" && <span className="stamp" style={{ marginLeft: 8 }}>sorted by outcome · exploratory</span>}</h2>
            <div className="seg">{["atlas", "mean"].map((o) => <button key={o} className={order === o ? "on" : ""} onClick={() => setOrder(o)}>{o === "atlas" ? "Atlas order" : "By mean"}</button>)}</div>
          </div>
          <RegionPairs regions={res.regions} subjects={res.subjects} chance={res.chance} order={order} hover={hover} setHover={setHover} />
          <p className="fig-caption">A dot per subject (hover one to find it in every row); the tick is the mean, the bar its 95% t-interval. These are {intervals} intervals with no adjustment for multiplicity, so one or two may clear chance by luck. Atlas order carries no selection.</p>
        </section>
      </div>
    );
  }

  const pz = res.peak_z;
  const slices = [-3, 0, 3, 6, 9, 12].map((d) => pz + d);
  const subjects = res.subjects.map((s) => ({ ...s, status: statusOf[s.id] || s.status }));
  // Open on the median subject (by LOSO accuracy), never a hand-picked one.
  const pool = subjects.filter((s) => ran.includes(s.id) && inN(s.id)).sort((p, q) => p.loso_roi - q.loso_roi);
  const median = pool[Math.floor((pool.length - 1) / 2)];
  const selected = subjects.find((s) => s.id === sel) || median;
  const tiles = [...a.cells].sort((p, q) => {
    if (tileOrder === "id") return p.subject.localeCompare(q.subject);
    const v = (c) => subjects.find((s) => s.id === c.subject)?.loso_roi ?? -1;
    return v(q) - v(p);
  });
  const prev = include ? maps.prevalence_all : maps.prevalence;
  const tops = Math.max(...Array.from(prev || []).filter((v) => v === v));
  const nRan = ran.filter((id) => inN(id)).length;
  const provisional = res.null.permutations < 1000;
  return (
    <div className="page plain">
      <h1>{a.name}</h1>
      <p className="lede">{res.metric} · {res.classes} classes, chance {res.chance}. <b>Inference is per subject:</b> each subject's searchlight map is tested against its own pooled permutation null ({res.null.adjust}). Views across subjects are descriptive; no group-level test is computed.</p>
      <Provenance a={a} res={res} used={n} excluded={excluded} />

      <section className="fig">
        <div className="fig-head">
          <h2>{view === "count" ? `How many of the ${nRan} subjects are significant at each centre` : `Mean accuracy of ${nRan} subjects`}</h2>
          <div className="seg">
            <button className={view === "count" ? "on" : ""} onClick={() => setView("count")}>Subjects significant</button>
            <button className={view === "mean" ? "on" : ""} onClick={() => setView("mean")}>Mean accuracy</button>
          </div>
        </div>
        {view === "count"
          ? <Slices template={vol.template} map={prev} dims={vol.dims} mode="count" top={tops} slices={slices} />
          : <Slices template={vol.template} map={maps.group_mean} dims={vol.dims} threshold={res.chance} top={0.42} slices={slices} />}
        <div className="subhead">
          <p className="fig-caption" style={{ margin: 0 }}>{view === "count"
            ? <>Descriptive: a count of independent per-subject decisions (own FDR, q &lt; .05), up to {tops} of {nRan}. There is no null for this count; a lone subject is often a false positive, so 1 is drawn faint.{provisional ? ` Provisional: ${res.null.permutations} permutations.` : ""}</>
            : <>Descriptive. Continuous from chance ({res.chance}) to 0.42; no threshold is implied.</>}</p>
          {view === "count" ? <CountKey bins={COUNT_BINS} /> : <AccuracyKey chance={res.chance} label="" />}
        </div>
      </section>

      <section className="fig">
        <div className="fig-head"><h2>sub-{selected.id} against its own null</h2><span className="origin">{sel ? "" : "median subject by LOSO accuracy · "}click any dot or tile · ← → to step</span></div>
        <SubjectNull subject={selected} edges={res.hist_edges} chance={res.chance} resolution={res.null.p_resolution} permutations={res.null.permutations} />
        <SubjectStrip subjects={subjects} chance={res.chance} include={include} hover={hover} setHover={setHover} select={setSel} selected={selected.id} />
        <div className="subhead">
          <p className="fig-caption" style={{ margin: 0 }}>Both panels share the accuracy axis. A subject's strip value is its accuracy in its leave-one-subject-out region: {res.group.roi_rule}, with <b>k = {res.group.roi_k}</b> <span className="origin">studio default, not yet in the analysis file or the lock</span>. No subject picks its own voxels. The interval describes the spread of these values; it is not a test.</p>
          <label className="toggle"><input id="include-stale" type="checkbox" checked={include} onChange={(e) => setInclude(e.target.checked)} /> include stale</label>
        </div>
      </section>

      <section className="fig">
        <div className="fig-head"><h2>Each subject at z {pz}: centres that survive its own FDR</h2>
          <div className="seg">{[["id", "By subject"], ["value", "By value"]].map(([k, l]) => <button key={k} className={tileOrder === k ? "on" : ""} onClick={() => setTileOrder(k)}>{l}</button>)}</div>
        </div>
        <div className="multiples">
          {tiles.map((c) => {
            const s = subjects.find((x) => x.id === c.subject);
            if (!s || !planes[c.subject]) return <div key={c.subject} className="tile missing">sub-{c.subject}<br />not run</div>;
            return (
              <div key={c.subject} className={`tile ${s.status === "stale" ? "stale" : ""} ${hover === c.subject ? "hl" : ""} ${selected.id === c.subject ? "sel" : ""}`}
                onMouseEnter={() => setHover(c.subject)} onMouseLeave={() => setHover(null)} onClick={() => setSel(c.subject)}>
                <PlaneSlice template={vol.template} plane={planes[c.subject]} dims={vol.dims} k={pz} threshold={s.fdr?.threshold ?? 9} mode="binary" />
                <div className="lab"><span><b>{c.subject}</b> · {s.loso_roi?.toFixed(3)}</span><span className={s.status === "stale" ? "stale-lab" : ""}>{s.status === "stale" ? "stale" : `${s.fdr?.surviving || 0} sig${provisional ? "*" : ""}`}</span></div>
              </div>
            );
          })}
        </div>
        <p className="fig-caption">Red marks the centres that survive the subject's own FDR; colour carries no accuracy, so tiles compare directly. Labels give the LOSO accuracy and the number of surviving centres{provisional ? ` (* provisional: ${res.null.permutations} permutations)` : ""}.</p>
      </section>
    </div>
  );
}
