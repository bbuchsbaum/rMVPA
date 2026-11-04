import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// One build serves every case (see src/api.js): rmvpa_studio() serves it from
// rmvpadsl's inst/app, rmvpa_export_studio() copies it next to the data, and
// `--mode example` adds the recorded example project from public/.
//
// The bundle is a single classic script (IIFE, dynamic imports inlined) loaded
// with `defer`: browsers refuse ES-module scripts on file:// pages, and an
// export has to open from a double-click.
const classicScript = () => ({
  name: "classic-script",
  enforce: "post",
  transformIndexHtml(html) {
    return html
      .replace(/<script type="module" crossorigin/g, "<script defer")
      .replace(/<link rel="stylesheet" crossorigin/g, '<link rel="stylesheet"');
  },
});

export default defineConfig(({ mode }) => ({
  plugins: [react(), classicScript()],
  base: "./",
  build: {
    outDir: mode === "example" ? "dist-example" : "../rmvpadsl/inst/app",
    emptyOutDir: true,
    copyPublicDir: mode === "example",
    modulePreload: false,
    chunkSizeWarningLimit: 1500,
    rollupOptions: { output: { format: "iife", inlineDynamicImports: true, entryFileNames: "assets/studio.js", assetFileNames: "assets/studio[extname]" } },
  },
  server: { proxy: { "/api": process.env.RMVPA_API || "http://127.0.0.1:7788" } },
}));
