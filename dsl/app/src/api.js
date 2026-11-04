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
export const api = {
  project: () => get("project", EX ? "./example/project.json" : "api/project"),
  template: () => get("template", EX ? "./example/template.json" : "api/template"),
  results: (name) => get(`results/${name}`, EX ? `./example/results-${name}.json` : `api/results?analysis=${q(name)}`),
  surface: (template) => get(`surface/${template}`, EX ? `./example/surface-${template}.json` : `api/surface?template=${q(template)}`),
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
