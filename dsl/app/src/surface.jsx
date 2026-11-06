// Parcel values on the cortical surface, drawn by surfview. Loaded lazily so
// three.js stays out of the rest of the app.
//
// Four small views (left lateral, left medial, right medial, right lateral),
// each its own NeuroSurfaceViewer. Per-vertex values come from neuroatlasjs's
// SurfaceAtlas.project; the medial wall (parcel 0) is left uncoloured.
import React, { useEffect, useRef } from "react";
import { NeuroSurfaceViewer, SurfaceGeometry, MultiLayerNeuroSurface, DataLayer, CurvatureLayer, OutlineLayer } from "surfview";
import { diverge } from "./brain.jsx";

const ZOOM = Number(new URLSearchParams(location.search).get("zoom")) || 255; // camera distance in mm for an fsaverage hemisphere
const VIEWS = [["left", "lateral"], ["left", "medial"], ["right", "medial"], ["right", "lateral"]];

function HemiView({ hemi, view, mesh, values, selected, span, deadzone, onPick, width, height }) {
  const box = useRef(null);
  const state = useRef(null);
  useEffect(() => {
    const viewer = new NeuroSurfaceViewer(box.current, width, height, { preset: "paper-light", showControls: false, backgroundColor: 0xf5f5f2 });
    const geometry = new SurfaceGeometry(mesh.positions, mesh.faces, hemi);
    const surface = new MultiLayerNeuroSurface(geometry, { baseColor: 0xd9d8d2 });
    surface.addLayer(new CurvatureLayer("sulc", mesh.sulc, {}));
    viewer.addSurface(surface, `${hemi}-${view}`);
    viewer.centerCamera();
    viewer.setViewpoint(view === "lateral" ? "lateral" : "medial");
    viewer.setZoom(ZOOM);
    viewer.on("surface:click", (hit) => { if (hit.vertexIndex != null) onPick?.(mesh.labels[hit.vertexIndex]); });
    state.current = { viewer, surface };
    return () => { viewer.dispose(); state.current = null; };
  }, [mesh, hemi, view, width, height]);

  useEffect(() => {
    const st = state.current;
    if (!st) return;
    for (const id of ["values", "selected"]) { try { st.surface.removeLayer(id); } catch { /* absent */ } }
    // Sample the same dead-zone ramp the slices use, so both views match.
    const ramp = Array.from({ length: 41 }, (_, i) => "#" + diverge(i / 20 - 1, deadzone).map((v) => v.toString(16).padStart(2, "0")).join(""));
    st.surface.addLayer(new DataLayer("values", values, null, ramp, { range: [-span, span] }));
    const present = (selected || []).filter((id) => mesh.labels.includes(id));
    if (present.length) {
      st.surface.addLayer(new OutlineLayer("selected", { roiLabels: mesh.labels, roiSubset: present, color: 0x17181b, width: 3, halo: false }));
    }
    st.viewer.requestRender?.();
  }, [values, (selected || []).join(), span, deadzone]);

  return <div ref={box} className="hemi-view" style={{ width, height }} />;
}

export default function CortexViews({ surface, projected, selected, span, deadzone, onPick, width = 250, height = 150 }) {
  return (
    <div className="cortex">
      {VIEWS.map(([h, v]) => (
        <figure key={h + v}>
          <HemiView hemi={h} view={v} mesh={surface[h]} values={projected[h]} selected={selected} span={span} deadzone={deadzone} onPick={onPick} width={width} height={height} />
          <figcaption>{h === "left" ? "Left" : "Right"} {v}</figcaption>
        </figure>
      ))}
    </div>
  );
}
