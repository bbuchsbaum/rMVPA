// Where the studio's data comes from, decided at run time so one build serves
// every case:
//   - a static export (rmvpa_export_studio) puts the data in window.__RMVPA_DATA__
//     through <script> tags, so the folder opens from file:// with no server;
//   - rmvpa_studio() serves the page with a session token in a <meta> tag and
//     answers /api/* only when that token comes back in a header;
//   - otherwise (vite dev, the hosted example) the recorded example project is
//     read from ./example/.
const DATA = typeof window !== "undefined" ? window.__RMVPA_DATA__ : undefined;
const TOKEN = typeof document !== "undefined" ? document.querySelector('meta[name="rmvpa-token"]')?.content : undefined;
export const source = DATA ? "export" : TOKEN ? "server" : "example";

const cache = new Map();
function get(key, url) {
  if (cache.has(key)) return cache.get(key);
  let p;
  if (DATA) {
    p = key in DATA ? Promise.resolve(DATA[key]) : Promise.reject(new Error(`${key}: not in this export`));
  } else {
    p = fetch(url, TOKEN ? { headers: { "X-Rmvpa-Token": TOKEN } } : undefined).then(async (r) => {
      if (!r.ok) throw new Error(`${url}: ${r.status} ${(await r.text()).slice(0, 200)}`);
      return r.json();
    });
  }
  cache.set(key, p);
  return p;
}

const q = encodeURIComponent;
const EX = source === "example";
const EXAMPLE_BACKGROUNDS = [{ id: "example", label: "MNI template (neuroim2, 4 mm)", status: "ready", default: true }];
export const api = {
  project: () => get("project", EX ? "./example/project.json" : "api/project"),
  // What can be drawn under the slices; an export carries only the one it used.
  backgrounds: () => (EX ? Promise.resolve(EXAMPLE_BACKGROUNDS) : get("backgrounds", "api/backgrounds")),
  template: (bg) => (EX || DATA ? get("template", EX ? "./example/template.json" : null)
    : get(`template/${bg || ""}`, `api/template${bg ? `?background=${q(bg)}` : ""}`)),
  results: (name) => get(`results/${name}`, EX ? `./example/results-${name}.json` : `api/results?analysis=${q(name)}`),
  surface: (id) => get(`surface/${id}`, EX ? `./example/surface-${id}.json` : `api/surface?id=${q(id)}`),
  // Asks the studio server to download (or compute) an asset it offered, then
  // forgets cached answers that may change because of it.
  canFetch: source === "server",
  fetch: async (body) => {
    if (source !== "server") throw new Error("Fetching needs the studio server: run `rmvpa studio`, or export with fetch = TRUE.");
    const r = await fetch("api/fetch", { method: "POST", headers: { "Content-Type": "application/json", "X-Rmvpa-Token": TOKEN }, body: JSON.stringify(body) });
    const out = await r.json().catch(() => ({}));
    if (!r.ok) throw new Error(out.error || `fetch failed (${r.status})`);
    for (const k of [...cache.keys()]) if (k === "backgrounds" || k.startsWith("results/")) cache.delete(k);
    return out;
  },
};

export function decodeF32(b64) {
  const bin = atob(b64);
  const buf = new ArrayBuffer(bin.length);
  const u8 = new Uint8Array(buf);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  return new Float32Array(buf);
}
export function decodeU8(b64) {
  const bin = atob(b64);
  const u8 = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  return u8;
}
