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
export const DIVERGING = [[40, 70, 170], [120, 150, 220], [232, 230, 224], [236, 150, 120], [212, 73, 44], [140, 35, 18]];
export function diverge(t) { // t in [-1, 1]
  t = (Math.max(-1, Math.min(1, t)) + 1) / 2 * (DIVERGING.length - 1);
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
      if (v === v && mode === "diverge") {
        c = diverge((v - threshold) / top);
        if (labelAt && opts.selected != null && labelAt(i, j) === opts.selected) {
          const L = (a, b) => labelAt(a, b);
          const edge = i === 0 || j === 0 || i === nx - 1 || j === ny - 1 ||
            L(i - 1, j) !== opts.selected || L(i + 1, j) !== opts.selected || L(i, j - 1) !== opts.selected || L(i, j + 1) !== opts.selected;
          if (edge) c = [23, 24, 27];
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
export function ParcelSlices({ template, values, labels, dims, center, span, slices, selected, onPick }) {
  return (
    <div className="slices lr" style={{ gridTemplateColumns: `repeat(${slices.length}, minmax(0, 1fr))` }}>
      {slices.map((k) => <ParcelSlice key={k} k={k} template={template} values={values} labels={labels} dims={dims}
        center={center} span={span} selected={selected} onPick={onPick} />)}
    </div>
  );
}

function ParcelSlice({ template, values, labels, dims, k, center, span, selected, onPick }) {
  const ref = useRef(null);
  const [nx, ny] = dims;
  const at = (i, j) => i + nx * (j + ny * k);
  useEffect(() => {
    paint(ref.current, nx, ny, (i, j) => template[at(i, j)], (i, j) => values[at(i, j)],
      { mode: "diverge", threshold: center, top: span, selected }, (i, j) => labels[at(i, j)]);
  }, [template, values, labels, k, center, span, selected]);
  const click = (e) => {
    const r = e.currentTarget.getBoundingClientRect();
    const i = Math.floor((e.clientX - r.left) / r.width * nx), j = ny - 1 - Math.floor((e.clientY - r.top) / r.height * ny);
    const id = labels[at(i, j)];
    if (id > 0) onPick?.(id);
  };
  return <div><canvas ref={ref} onClick={click} style={{ cursor: "pointer" }} /><div className="z">z {k}</div></div>;
}
