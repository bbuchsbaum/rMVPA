import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// The built app ships inside rmvpadsl (inst/app) and is served by the R
// package. `--mode example` builds a static copy that reads the recorded
// example project from ./example/ instead of the R API.
export default defineConfig(({ mode }) => ({
  plugins: [react()],
  base: "./",
  define: { __EXAMPLE__: JSON.stringify(mode === "example") },
  build: {
    outDir: mode === "example" ? "dist-example" : "../rmvpadsl/inst/app",
    emptyOutDir: true,
    copyPublicDir: mode === "example",
  },
  server: { proxy: { "/api": process.env.RMVPA_API || "http://127.0.0.1:7788" } },
}));
