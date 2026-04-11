(function () {
  var FAMILY_CLASSES = ["red", "lapis", "ochre", "teal", "green", "violet"];
  var PRESET_CLASSES = ["study", "structural", "adobe", "midnight"];
  var STYLE_CLASSES = ["minimal", "assertive"];

  function classes(prefix, values) {
    return values.map(function (v) { return prefix + v; });
  }

  function applySingleClass(el, prefix, value, values) {
    classes(prefix, values).forEach(function (c) { el.classList.remove(c); });
    if (value) el.classList.add(prefix + value);
  }

  function inferBootstrapTheme() {
    if (!document.body) return "light";
    return document.body.classList.contains("preset-midnight") ? "dark" : "light";
  }

  function syncBootstrapTheme() {
    var theme = inferBootstrapTheme();
    document.documentElement.setAttribute("data-bs-theme", theme);
    document.body.setAttribute("data-bs-theme", theme);
    var nav = document.querySelector("nav.navbar");
    if (nav) nav.setAttribute("data-bs-theme", theme);
  }

  function hashSeed(seed) {
    var h = 2166136261;
    for (var i = 0; i < seed.length; i++) {
      h ^= seed.charCodeAt(i);
      h += (h << 1) + (h << 4) + (h << 7) + (h << 8) + (h << 24);
    }
    return h >>> 0;
  }

  function seeded(seed) {
    var state = hashSeed(seed || "albersdown");
    return function () {
      state = (1664525 * state + 1013904223) >>> 0;
      return state / 4294967296;
    };
  }

  function readPalette(el) {
    var source = el || document.body;
    var style = getComputedStyle(source);
    var colors = ["--A900", "--A700", "--A500", "--A300"].map(function (k) {
      var value = style.getPropertyValue(k).trim();
      return value || "#666";
    });
    return colors;
  }

  function clearChildren(el) {
    while (el.firstChild) el.removeChild(el.firstChild);
  }

  function svgEl(name, attrs) {
    var el = document.createElementNS("http://www.w3.org/2000/svg", name);
    Object.keys(attrs || {}).forEach(function (key) {
      el.setAttribute(key, String(attrs[key]));
    });
    return el;
  }

  function renderComposition(container, idx) {
    var seed = container.getAttribute("data-seed") || (window.location.pathname + ":" + idx);
    var density = Number(container.getAttribute("data-density") || "6");
    var rand = seeded(seed);
    var palette = readPalette(container);
    var w = 1600;
    var h = 500;

    clearChildren(container);
    var svg = svgEl("svg", {
      viewBox: "0 0 " + w + " " + h,
      preserveAspectRatio: "xMidYMid slice",
      role: "img",
      "aria-label": "Albers composition"
    });

    svg.appendChild(svgEl("rect", { x: 0, y: 0, width: w, height: h, fill: "var(--surface)" }));

    var blockCount = Math.max(4, Math.min(16, density * 2));
    for (var i = 0; i < blockCount; i++) {
      var side = 90 + rand() * 260;
      var x = rand() * (w - side);
      var y = rand() * (h - side);
      var tone = palette[Math.floor(rand() * palette.length)];
      svg.appendChild(svgEl("rect", {
        x: x.toFixed(2),
        y: y.toFixed(2),
        width: side.toFixed(2),
        height: side.toFixed(2),
        fill: tone,
        opacity: (0.17 + rand() * 0.72).toFixed(3)
      }));
    }

    for (var j = 0; j < 3; j++) {
      var cx = 80 + rand() * (w - 160);
      var cy = 60 + rand() * (h - 120);
      var radius = 26 + rand() * 66;
      var color = palette[Math.floor(rand() * palette.length)];
      svg.appendChild(svgEl("circle", {
        cx: cx.toFixed(2),
        cy: cy.toFixed(2),
        r: radius.toFixed(2),
        fill: color,
        opacity: (0.16 + rand() * 0.44).toFixed(3)
      }));
    }

    svg.appendChild(svgEl("rect", {
      x: 0,
      y: 0,
      width: w,
      height: h,
      fill: "none",
      stroke: "var(--border)",
      "stroke-width": "4"
    }));

    container.appendChild(svg);
  }

  function renderAllCompositions() {
    Array.from(document.querySelectorAll(".albers-composition")).forEach(function (el, idx) {
      renderComposition(el, idx);
    });
  }

  function initCopyButtons() {
    var hasPkgdownCopy =
      document.querySelector(".btn-copy-ex") ||
      document.querySelector("[data-clipboard-copy]") ||
      document.querySelector("div.sourceCode.hasCopyButton");

    if (hasPkgdownCopy) return;

    document.querySelectorAll("pre > code").forEach(function (code) {
      var pre = code.parentElement;
      if (!pre || pre.querySelector("button.copy-code")) return;
      if (!code.textContent || !code.textContent.trim()) return;
      var btn = document.createElement("button");
      btn.className = "copy-code";
      btn.type = "button";
      btn.setAttribute("aria-label", "Copy code");
      btn.textContent = "Copy";
      btn.addEventListener("click", async function () {
        try {
          await navigator.clipboard.writeText(code.innerText);
          btn.textContent = "Copied";
        } catch (e) {
          try {
            var r = document.createRange();
            r.selectNodeContents(code);
            var sel = window.getSelection();
            sel.removeAllRanges();
            sel.addRange(r);
            document.execCommand("copy");
            sel.removeAllRanges();
            btn.textContent = "Copied";
          } catch (e2) {
            btn.textContent = "Oops";
          }
        } finally {
          setTimeout(function () { btn.textContent = "Copy"; }, 1400);
        }
      });
      pre.appendChild(btn);
    });
  }

  function initAnchors() {
    Array.from(document.querySelectorAll("h2, h3")).forEach(function (h) {
      if (!h.id || h.querySelector("a.anchor")) return;
      var a = document.createElement("a");
      a.href = "#" + h.id;
      a.className = "anchor";
      a.textContent = "—";
      a.setAttribute("aria-label", "Link to this section");
      a.setAttribute("title", "Link to this section");
      h.appendChild(a);
    });
  }

  function initDefaults() {
    if (
      document.body &&
      !document.body.classList.contains("style-minimal") &&
      !document.body.classList.contains("style-assertive")
    ) {
      document.body.classList.add("style-minimal");
    }
  }

  function setSearchParam(key, value) {
    try {
      var url = new URL(window.location.href);
      if (value === null || value === "") url.searchParams.delete(key);
      else url.searchParams.set(key, value);
      window.history.replaceState({}, "", url.toString());
    } catch (e) {
      // no-op
    }
  }

  function updateLabSummary(root, state) {
    var target = root.querySelector("[data-albers-lab-summary]");
    if (!target) return;
    target.textContent = [
      "family=" + state.family,
      "preset=" + state.preset,
      "style=" + state.style,
      "width=" + state.width + "ch"
    ].join(" | ");
  }

  function applyLabState(root, state) {
    applySingleClass(document.body, "palette-", state.family, FAMILY_CLASSES);
    applySingleClass(document.body, "preset-", state.preset, PRESET_CLASSES);

    classes("style-", STYLE_CLASSES).forEach(function (c) {
      document.body.classList.remove(c);
    });
    if (state.style && state.style !== "balanced") {
      document.body.classList.add("style-" + state.style);
    }

    document.documentElement.style.setProperty("--content", state.width + "ch");

    syncBootstrapTheme();

    setSearchParam("family", state.family);
    setSearchParam("preset", state.preset);
    setSearchParam("style", state.style);
    setSearchParam("width", String(state.width));

    updateLabSummary(root, state);
    renderAllCompositions();
  }

  function watchThemeClasses() {
    if (!document.body || typeof MutationObserver === "undefined") return;
    var observer = new MutationObserver(function (muts) {
      for (var i = 0; i < muts.length; i++) {
        if (muts[i].attributeName === "class") {
          syncBootstrapTheme();
          renderAllCompositions();
          break;
        }
      }
    });
    observer.observe(document.body, { attributes: true, attributeFilter: ["class"] });
  }

  function initThemeLab() {
    var root = document.querySelector("[data-albers-lab]");
    if (!root) return;

    var controls = {
      family: root.querySelector("[data-albers-control='family']"),
      preset: root.querySelector("[data-albers-control='preset']"),
      style: root.querySelector("[data-albers-control='style']"),
      width: root.querySelector("[data-albers-control='width']")
    };

    if (!controls.family || !controls.preset || !controls.style || !controls.width) return;

    var params = new URLSearchParams(window.location.search);
    var initial = {
      family: params.get("family") || controls.family.value || "red",
      preset: params.get("preset") || controls.preset.value || "study",
      style: params.get("style") || controls.style.value || "minimal",
      width: Number(params.get("width") || controls.width.value || 78)
    };

    controls.family.value = initial.family;
    controls.preset.value = initial.preset;
    controls.style.value = initial.style;
    controls.width.value = String(initial.width);

    applyLabState(root, initial);

    Object.keys(controls).forEach(function (key) {
      controls[key].addEventListener("change", function () {
        applyLabState(root, {
          family: controls.family.value,
          preset: controls.preset.value,
          style: controls.style.value,
          width: Number(controls.width.value || "78")
        });
      });
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    initDefaults();
    initCopyButtons();
    initAnchors();
    initThemeLab();
    watchThemeClasses();
    syncBootstrapTheme();
    renderAllCompositions();
  });
})();
