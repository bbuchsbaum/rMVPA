// The studio talks to rmvpadsl's local server (/api/*). The example build
// reads a recorded project from ./example/ instead.
/* global __EXAMPLE__ */
const EXAMPLE = typeof __EXAMPLE__ !== "undefined" && __EXAMPLE__;

const cache = new Map();
function get(url) {
  if (!cache.has(url)) cache.set(url, fetch(url).then((r) => {
    if (!r.ok) throw new Error(`${url}: ${r.status}`);
    return r.json();
  }));
  return cache.get(url);
}

export const isExample = EXAMPLE;
export const api = {
  project: () => get(EXAMPLE ? "./example/project.json" : "/api/project"),
  template: () => get(EXAMPLE ? "./example/template.json" : "/api/template"),
  results: (name) => get(EXAMPLE ? `./example/results-${name}.json` : `/api/results?analysis=${encodeURIComponent(name)}`),
  surface: (template) => get(EXAMPLE ? `./example/surface-${template}.json` : `/api/surface?template=${encodeURIComponent(template)}`),
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
