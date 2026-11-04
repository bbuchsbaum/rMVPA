import React, { useEffect, useRef } from "react";

// Axial slices of a result map over the template, drawn on canvas.
// Arrays are x-fastest (R order): i + nx * (j + ny * k).
// Rows run anterior (top) to posterior; columns run left to right (neurological).
// Values at or above `threshold` take the hot ramp; values between `floor`
// and `threshold` are tinted faintly, so sub-threshold signal stays visible.

const HOT = [[122, 29, 14], [212, 73, 44], [243, 177, 59], [255, 241, 194]];
export function hot(t) {
  t = Math.max(0, Math.min(1, t)) * (HOT.length - 1);
  const i = Math.min(HOT.length - 2, Math.floor(t)), f = t - i;
  return HOT[i].map((c, k) => Math.round(c + f * (HOT[i + 1][k] - c)));
}

// mode "count": integer counts (e.g. subjects significant); 0 shows the template.
function paint(cv, nx, ny, tplAt, valAt, { threshold, top, floor, mode }) {
  cv.width = nx; cv.height = ny;
  const ctx = cv.getContext("2d");
  const img = ctx.createImageData(nx, ny);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const p = ((ny - 1 - j) * nx + i) * 4;
      const g = tplAt(i, j);
      let c = [g, g, g];
      const v = valAt(i, j);
      if (v === v && mode === "count") {
        if (v >= 2) c = hot((v - 2) / Math.max(1, top - 2));
        else if (v >= 1) c = [g + 0.3 * (212 - g), g + 0.3 * (73 - g), g + 0.3 * (44 - g)];
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
  const ref = useRef(null);
  useEffect(() => {
    const [nx, ny] = dims;
    paint(ref.current, nx, ny, (i, j) => template[i + nx * (j + ny * k)], (i, j) => plane[i + nx * j], opts);
  }, [template, plane, dims, k, opts.threshold, opts.top, opts.floor]);
  return <canvas ref={ref} />;
}
