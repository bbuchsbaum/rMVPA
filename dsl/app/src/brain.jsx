import React, { useEffect, useRef } from "react";

// Axial slices of a result map over the template, drawn on canvas.
// Arrays are x-fastest (R order): i + nx * (j + ny * k).
// Rows run anterior (top) to posterior; columns run left to right (neurological).
// Values at or above `threshold` take the hot ramp; values between `floor`
// and `threshold` are tinted faintly, so sub-threshold signal stays visible.

// Subjects-significant counts, binned so each step is clearly distinct.
export const COUNT_BINS = [
  { lo: 1, hi: 1, label: "1", rgb: [214, 190, 182] },
  { lo: 2, hi: 3, label: "2–3", rgb: [150, 40, 22] },
  { lo: 4, hi: 6, label: "4–6", rgb: [224, 92, 44] },
  { lo: 7, hi: 99, label: "7+", rgb: [250, 214, 120] },
];

// Diverging ramp around chance: blue below, neutral at chance, warm above.
// Odd number of stops so the centre (chance) is exactly the neutral one.
export const DIVERGING = [[40, 70, 170], [95, 125, 210], [170, 185, 230], [236, 234, 229], [240, 165, 135], [212, 73, 44], [140, 35, 18]];
// A small neutral band around chance (dz, as a fraction of the span) keeps
// near-chance parcels from reading as "above chance".
export function diverge(t, dz = 0) { // t in [-1, 1]
  t = Math.max(-1, Math.min(1, t));
  t = Math.sign(t) * Math.max(0, (Math.abs(t) - dz) / (1 - dz));
  t = (t + 1) / 2 * (DIVERGING.length - 1);
  const i = Math.min(DIVERGING.length - 2, Math.floor(t)), f = t - i;
  return DIVERGING[i].map((c, k) => Math.round(c + f * (DIVERGING[i + 1][k] - c)));
}

const HOT = [[122, 29, 14], [212, 73, 44], [243, 177, 59], [255, 241, 194]];
export function hot(t) {
  t = Math.max(0, Math.min(1, t)) * (HOT.length - 1);
  const i = Math.min(HOT.length - 2, Math.floor(t)), f = t - i;
  return HOT[i].map((c, k) => Math.round(c + f * (HOT[i + 1][k] - c)));
}

// mode "count": integer counts (e.g. subjects significant); 0 shows the template.
function paint(cv, nx, ny, tplAt, valAt, opts, labelAt) {
  const { threshold, top, floor, mode } = opts;
  cv.width = nx; cv.height = ny;
  const ctx = cv.getContext("2d");
  const img = ctx.createImageData(nx, ny);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const p = ((ny - 1 - j) * nx + i) * 4;
      const g = tplAt(i, j);
      let c = [g, g, g];
      const v = valAt(i, j);
      if (mode === "diverge") {
        // Paper ground, faint anatomy, parcel colour, ink outline with a paper halo.
        const paperC = [245, 245, 242];
        c = g < 8 ? paperC : paperC.map((p) => p - (p - g) * 0.35);
        if (v === v) c = diverge((v - threshold) / top, opts.deadzone || 0);
        const sel = opts.selected || [];
        if (labelAt && sel.length) {
          const me = labelAt(i, j);
          const nb = [[i - 1, j], [i + 1, j], [i, j - 1], [i, j + 1]].map(([a, b]) => (a < 0 || b < 0 || a >= nx || b >= ny ? 0 : labelAt(a, b)));
          if (sel.includes(me) && nb.some((n) => n !== me)) c = [23, 24, 27];
          else if (!sel.includes(me) && nb.some((n) => sel.includes(n))) c = paperC;
        }
      } else if (v === v && mode === "binary") {
        if (v >= threshold) c = [212, 73, 44];
      } else if (v === v && mode === "count") {
        const bin = COUNT_BINS.find((b) => v >= b.lo && v <= b.hi);
        if (bin) c = bin.rgb;
      } else if (v === v && threshold != null) {
        if (v >= threshold) c = hot((v - threshold) / Math.max(1e-9, top - threshold));
        else if (floor != null && v > floor) {
          const a = 0.45 * Math.min(1, (v - floor) / Math.max(1e-9, threshold - floor));
          c = [g + a * (212 - g), g + a * (73 - g), g + a * (44 - g)];
        }
      }
      img.data[p] = c[0]; img.data[p + 1] = c[1]; img.data[p + 2] = c[2]; img.data[p + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
}

export function Slices({ template, map, dims, threshold, top, floor, slices, mode }) {
  return (
    <div className="slices lr" style={{ gridTemplateColumns: `repeat(${slices.length}, minmax(0, 1fr))` }}>
      {slices.map((k) => (
        <div key={k}>
          <VolSlice template={template} map={map} dims={dims} k={k} threshold={threshold} top={top} floor={floor} mode={mode} />
          <div className="z">z {k}</div>
        </div>
      ))}
    </div>
  );
}

function VolSlice({ template, map, dims, k, ...opts }) {
  const ref = useRef(null);
  useEffect(() => {
    const [nx, ny] = dims;
    paint(ref.current, nx, ny, (i, j) => template[i + nx * (j + ny * k)], (i, j) => (map ? map[i + nx * (j + ny * k)] : NaN), opts);
  }, [template, map, dims, k, opts.threshold, opts.top, opts.floor, opts.mode]);
  return <canvas ref={ref} />;
}

// One 2D slice (x-fastest) over the template slice k.
export function PlaneSlice({ template, plane, dims, k, ...opts }) {
  // opts: threshold, top, floor, mode
  const ref = useRef(null);
  useEffect(() => {
    const [nx, ny] = dims;
    paint(ref.current, nx, ny, (i, j) => template[i + nx * (j + ny * k)], (i, j) => plane[i + nx * j], opts);
  }, [template, plane, dims, k, opts.threshold, opts.top, opts.floor, opts.mode]);
  return <canvas ref={ref} />;
}

// Parcel values painted into the atlas on axial slices. `values` is a
// per-voxel Float64Array (NaN outside parcels) from VolumeAtlas.project;
// clicking a voxel reports its parcel id.
export function ParcelSlices({ template, values, labels, dims, center, span, deadzone, slices, selected, names, onPick }) {
  return (
    <div className="pslices" style={{ gridTemplateColumns: `repeat(${slices.length}, minmax(0, 1fr))` }}>
      {slices.map((k, n) => <ParcelSlice key={k} k={k} first={n === 0} template={template} values={values} labels={labels} dims={dims}
        center={center} span={span} deadzone={deadzone} selected={selected} names={names} onPick={onPick} />)}
    </div>
  );
}

function ParcelSlice({ template, values, labels, dims, k, first, center, span, deadzone, selected, names, onPick }) {
  const ref = useRef(null);
  const [nx, ny] = dims;
  const at = (i, j) => i + nx * (j + ny * k);
  useEffect(() => {
    paint(ref.current, nx, ny, (i, j) => template[at(i, j)], (i, j) => values[at(i, j)],
      { mode: "diverge", threshold: center, top: span, deadzone, selected }, (i, j) => labels[at(i, j)]);
  }, [template, values, labels, k, center, span, deadzone, (selected || []).join()]);
  // Name the selected parcel in place, at its centroid on this slice.
  const tags = (selected || []).map((id) => {
    let sx = 0, sy = 0, n = 0;
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) if (labels[at(i, j)] === id) { sx += i; sy += j; n++; }
    return n >= 4 ? { id, x: (sx / n + 0.5) / nx * 100, y: (1 - (sy / n + 0.5) / ny) * 100 } : null;
  }).filter(Boolean);
  const click = (e) => {
    const r = e.currentTarget.getBoundingClientRect();
    const i = Math.floor((e.clientX - r.left) / r.width * nx), j = ny - 1 - Math.floor((e.clientY - r.top) / r.height * ny);
    const id = labels[at(i, j)];
    if (id > 0) onPick?.(id);
  };
  return (
    <div className="pslice">
      <canvas ref={ref} onClick={click} style={{ cursor: "pointer" }} />
      {first && <><span className="lr-tag l">L</span><span className="lr-tag r">R</span></>}
      {tags.map((t) => <span key={t.id} className="ptag" style={{ left: `${t.x}%`, top: `${t.y}%` }}>{names?.[t.id]}</span>)}
      <div className="z">z {k}</div>
    </div>
  );
}
