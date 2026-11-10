// Affordances for images the studio finds or fetches, so nobody types a path:
// the background under the slices, and the cortical surface of a regional
// analysis. The server says what exists, what can be fetched and from where.
import React, { useEffect, useState } from "react";
import { api } from "./api.js";

const BG_KEY = "rmvpa.background";
export function useBackground() {
  const [bg, setBg] = useState(() => { try { return localStorage.getItem(BG_KEY) || null; } catch { return null; } });
  const choose = (id) => { setBg(id); try { localStorage.setItem(BG_KEY, id); } catch { /* private mode */ } };
  return [bg, choose];
}

// "Under the slices: [select]". A fetchable entry downloads on pick and then
// becomes the choice; a failure is said in place, the old choice stays.
export function BackgroundPicker({ value, onChange }) {
  const [list, setList] = useState(null);
  const [busy, setBusy] = useState(null);
  const [err, setErr] = useState(null);
  useEffect(() => { api.backgrounds().then(setList).catch(() => setList([])); }, []);
  if (!list || list.length < 2) {
    const only = list?.[0];
    return only ? <span className="bg-pick"><span className="bg-label">Background</span> <span className="bg-only">{only.label}</span></span> : null;
  }
  const current = list.find((b) => b.id === value && b.status === "ready") || list.find((b) => b.default) || list[0];
  const pick = async (id) => {
    const b = list.find((x) => x.id === id);
    setErr(null);
    if (b.status === "fetch") {
      setBusy(b.label);
      try { const out = await api.fetch({ kind: "background", id }); setList(out.backgrounds); onChange(id); }
      catch (e) { setErr(e.message); }
      finally { setBusy(null); }
    } else onChange(id);
  };
  return (
    <span className="bg-pick">
      <label className="bg-label" htmlFor="bg-select">Background</label>
      <select id="bg-select" className="select small" title={current.source} value={current.id} disabled={!!busy} onChange={(e) => pick(e.target.value)}>
        {list.map((b) => <option key={b.id} value={b.id}>{b.label}{b.status === "fetch" ? ` — fetch from ${b.source} (${b.size})` : ""}</option>)}
      </select>
      {busy && <span className="bg-note">Fetching {busy}…</span>}
      {err && <span className="bg-note bad" role="alert">{err}</span>}
    </span>
  );
}

// Shown in place of the cortex until its surface exists: what will be
// fetched, from where, and whether the mapping is exact.
export function SurfaceOffer({ option, analysis, onReady }) {
  const [busy, setBusy] = useState(false);
  const [err, setErr] = useState(null);
  if (!option || option.status === "unavailable") {
    return <div className="offer"><b>No cortical view for this atlas.</b> {option?.reason || ""} The slices show the atlas the analysis used.</div>;
  }
  const go = async () => {
    setBusy(true); setErr(null);
    try { await api.fetch({ kind: "surface", analysis }); onReady(); }
    catch (e) { setErr(e.message); setBusy(false); }
  };
  return (
    <div className="offer">
      <div>
        <b>{option.label}</b>{option.exact ? "" : <span className="tag-approx">approximate</span>}
        <p>{option.source}.</p>
        {!api.canFetch && <p className="bg-note">This is a static export: run <code>rmvpa studio</code> to fetch, or export with <code>fetch = TRUE</code>.</p>}
        {err && <p className="bg-note bad" role="alert">{err}</p>}
      </div>
      {api.canFetch && <button className="button primary" onClick={go} disabled={busy}>{busy ? "Fetching…" : `Fetch and draw (${option.size})`}</button>}
    </div>
  );
}
