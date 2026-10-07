(function () {
  // Loaded once per page (a site may link albers.js from the template and
  // also carry a copy in pkgdown/extra.js).
  if (window.__albersdown) return;
  window.__albersdown = true;

  // pkgdown links the site's extra.css before the template's albers.css;
  // move it after, so extra.css stays the override layer it is meant to be.
  (function () {
    var themed = document.querySelector('link[href$="albers.css"]');
    var extra = document.querySelector('link[href$="extra.css"]');
    if (themed && extra && (themed.compareDocumentPosition(extra) & Node.DOCUMENT_POSITION_PRECEDING)) {
      themed.parentNode.insertBefore(extra, themed.nextSibling);
    }
  })();

  var FAMILY_CLASSES = ["red", "lapis", "ochre", "teal", "green", "violet"];
  var PRESET_CLASSES = ["homage", "interaction", "study", "structural", "adobe", "midnight"];
  var STYLE_CLASSES = ["minimal", "balanced", "assertive"];
  var THEME_KEY = "albersdown-theme";
  var root = document.documentElement;

  function classes(prefix, values) {
    return values.map(function (v) { return prefix + v; });
  }

  function applySingleClass(el, prefix, value, values) {
    classes(prefix, values).forEach(function (c) { el.classList.remove(c); });
    if (value) el.classList.add(prefix + value);
  }

  /* ---------------------------------------------------------------------
   * Theme (light / dark). Runs immediately in <head> so the first paint is
   * already correct. Priority: stored choice > OS preference. The midnight
   * preset is always dark. pkgdown's own light-switch, when present, owns
   * [data-bs-theme]; we mirror it.
   * ------------------------------------------------------------------- */
  function storedTheme() {
    try { return localStorage.getItem(THEME_KEY); } catch (e) { return null; }
  }

  function storeTheme(value) {
    try {
      if (value) localStorage.setItem(THEME_KEY, value);
      else localStorage.removeItem(THEME_KEY);
    } catch (e) { /* storage unavailable: choice lasts for this page view */ }
  }

  function systemDark() {
    return !!(window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches);
  }

  function hasPkgdownLightswitch() {
    return !!document.querySelector("#dropdown-lightswitch, [data-bs-theme-value]");
  }

  var sessionChoice = null;

  function resolveTheme() {
    if (document.body && document.body.classList.contains("preset-midnight")) return "dark";
    var choice = sessionChoice || storedTheme();
    if (choice === "light" || choice === "dark") return choice;
    return systemDark() ? "dark" : "light";
  }

  function applyTheme() {
    if (document.body && hasPkgdownLightswitch()) {
      // pkgdown's light-switch is in charge; follow whatever it decided.
      var bs = root.getAttribute("data-bs-theme");
      root.setAttribute("data-albers-theme", bs === "dark" ? "dark" : "light");
      return;
    }
    var theme = resolveTheme();
    root.setAttribute("data-albers-theme", theme);
    root.setAttribute("data-bs-theme", theme);
    if (document.body) {
      document.body.setAttribute("data-bs-theme", theme);
      var nav = document.querySelector("nav.navbar");
      if (nav) nav.setAttribute("data-bs-theme", theme);
    }
    var btn = document.querySelector(".albers-theme-toggle");
    if (btn) {
      var mode = themeMode();
      var sys = systemDark() ? "dark" : "light";
      var nextLabel = mode === "auto" ? (theme === "dark" ? "light" : "dark")
        : (mode === sys ? "follow your system" : (mode === "dark" ? "light" : "dark"));
      var label = (mode === "auto" ? "Theme: follows your system (" + theme + ")" : "Theme: " + mode) +
        ". Switch to " + nextLabel;
      btn.setAttribute("aria-label", label);
      btn.setAttribute("title", label);
      btn.setAttribute("data-mode", mode);
      btn.innerHTML = TOGGLE_ICONS[mode];
    }
  }

  function themeMode() {
    var choice = sessionChoice || storedTheme();
    return choice === "light" || choice === "dark" ? choice : "auto";
  }

  // An Albers square in three states: empty core (light), filled core
  // (dark), and a core split down the middle (follow the system).
  var FRAME = '<rect x="1.2" y="1.2" width="13.6" height="13.6" fill="none" stroke="currentColor" stroke-width="1.4"/>';
  var TOGGLE_ICONS = {
    light: '<svg viewBox="0 0 16 16" aria-hidden="true" focusable="false">' + FRAME +
      '<rect x="5.1" y="6.9" width="5.8" height="5.8" fill="none" stroke="currentColor" stroke-width="1.4"/></svg>',
    dark: '<svg viewBox="0 0 16 16" aria-hidden="true" focusable="false">' + FRAME +
      '<rect x="4.4" y="6.2" width="7.2" height="7.2" fill="currentColor"/></svg>',
    auto: '<svg viewBox="0 0 16 16" aria-hidden="true" focusable="false">' + FRAME +
      '<rect x="5.1" y="6.9" width="5.8" height="5.8" fill="none" stroke="currentColor" stroke-width="1.4"/>' +
      '<rect x="8" y="6.2" width="3.6" height="7.2" fill="currentColor"/></svg>'
  };

  /* Site defaults (pkgdown/extra.js sets window.albersdownDefaults). Applied as
     soon as <body> is created -- before its content paints -- and only where
     no class of that kind is present yet. */
  function hasAny(el, prefix, values) {
    return values.some(function (v) { return el.classList.contains(prefix + v); });
  }

  function applySiteDefaults() {
    var d = window.albersdownDefaults;
    var b = document.body;
    if (!d || !b) return;
    if (d.family && !hasAny(b, "palette-", FAMILY_CLASSES)) b.classList.add("palette-" + d.family);
    if (d.preset && !hasAny(b, "preset-", PRESET_CLASSES)) b.classList.add("preset-" + d.preset);
    if (d.style && !hasAny(b, "style-", STYLE_CLASSES)) b.classList.add("style-" + d.style);
  }

  function whenBody(fn) {
    if (document.body) { fn(); return; }
    if (typeof MutationObserver === "undefined") {
      document.addEventListener("DOMContentLoaded", fn);
      return;
    }
    var mo = new MutationObserver(function () {
      if (document.body) { mo.disconnect(); fn(); }
    });
    mo.observe(document.documentElement, { childList: true });
  }

  whenBody(function () { applySiteDefaults(); applyTheme(); });

  /* Pages get their direction/family classes from an inline script in the
     body (vignettes, pkgdown articles), and albers.js then gathers the title
     block and figures. Hold the first paint until both have happened, so the
     page never shows the wrong theme or reflows (failsafe: 1.2s). */
  // A long page returns the reader to their place itself on a reload or
  // Back/Forward (Chrome's restoration of the saved pixel offset misfired
  // there, by up to 83,000 px): the line or element at the top of the window
  // is saved when the page is left, and restored once the code before it is
  // prepared.
  // (The browser still restores its own offset -- same-document entries and
  // pages without storage keep working -- and albers.js corrects it from
  // the saved place once the code before it is prepared. The key is per
  // history entry, so Back to an earlier visit gets that visit's place.)
  var PLACE_KEY = null, savedPlace = null;
  try {
    // (another script's state that is not a plain object is left alone;
    // the place then lives in session storage only)
    var st0 = history.state;
    var plainState = st0 === null || (typeof st0 === "object" && !Array.isArray(st0) && Object.getPrototypeOf(st0) === Object.prototype);
    var st = plainState ? st0 : null;
    var entry = st && st.albersEntry;
    if (!entry && plainState) {
      entry = String(Date.now()) + Math.random().toString(36).slice(2, 7);
      history.replaceState(Object.assign({}, st || {}, { albersEntry: entry }), "");
    }
    entry = entry || "page";
    PLACE_KEY = "albers-place:" + location.pathname + ":" + entry;
    var navEntry = performance.getEntriesByType ? performance.getEntriesByType("navigation")[0] : null;
    if (navEntry && navEntry.type !== "navigate") {
      // (the history entry's own copy first: it works where storage is
      // blocked)
      savedPlace = (st && st.albersPlace) || null;
      if (!savedPlace) { try { savedPlace = JSON.parse(sessionStorage.getItem(PLACE_KEY) || "null"); } catch (e2) {} }
    }
  } catch (e) { savedPlace = null; }

  var holding = true;
  if (holding) root.classList.add("albers-loading");
  function reveal() {
    if (!holding) return;
    holding = false;
    root.classList.remove("albers-loading");
  }
  setTimeout(reveal, 1200);

  applyTheme();

  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    var onChange = function () { if (!sessionChoice && !storedTheme()) applyTheme(); };
    if (mq.addEventListener) mq.addEventListener("change", onChange);
    else if (mq.addListener) mq.addListener(onChange);
  }

  function initThemeToggle(isVignette) {
    if (hasPkgdownLightswitch() || document.querySelector(".albers-theme-toggle")) return;
    if (document.body.classList.contains("preset-midnight")) return;
    var btn = document.createElement("button");
    btn.type = "button";
    btn.className = "albers-theme-toggle";
    btn.addEventListener("click", function () {
      if (document.body.classList.contains("preset-midnight")) return;
      var mode = themeMode();
      var now = resolveTheme();
      var next = mode === "auto" ? (now === "dark" ? "light" : "dark")
        : (mode === (systemDark() ? "dark" : "light") ? "auto" : (mode === "dark" ? "light" : "dark"));
      sessionChoice = next === "auto" ? null : next;
      storeTheme(sessionChoice);
      applyTheme();
      renderAllCompositions();
    });
    if (isVignette) {
      var byline = document.querySelector(".albers-titleblock .albers-byline");
      var tb = document.querySelector(".albers-titleblock");
      if (!byline && tb) {
        // no author/date: an empty byline row still holds the control
        byline = document.createElement("div");
        byline.className = "albers-byline";
        var plate = tb.querySelector(".albers-plate-mark");
        tb.insertBefore(byline, plate);
      }
      if (byline) {
        btn.classList.add("is-inline");
        byline.appendChild(btn);
      } else {
        document.body.appendChild(btn);
      }
    } else {
      // beside the menu button, so it stays visible when the menu collapses
      var toggler = document.querySelector(".navbar .navbar-toggler");
      if (toggler) toggler.parentNode.insertBefore(btn, toggler);
      else {
        var host = document.querySelector("#navbar") || document.querySelector(".navbar .container");
        if (!host) return;
        host.appendChild(btn);
      }
    }
    applyTheme();
  }

  /* ---------------------------------------------------------------------
   * Page type
   * ------------------------------------------------------------------- */
  function isPkgdown() {
    return !!document.querySelector('script[src*="pkgdown.js"], meta[name="generator"][content*="pkgdown"]');
  }

  function isVignetteDom() {
    return !!document.querySelector("body > h1.title, body > #TOC");
  }

  /* ---------------------------------------------------------------------
   * Generative composition motif
   * ------------------------------------------------------------------- */
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
    return ["--A900", "--A700", "--A500", "--A300"].map(function (k) {
      var value = style.getPropertyValue(k).trim();
      return value || "#666";
    });
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
      x: 0, y: 0, width: w, height: h,
      fill: "none", stroke: "var(--border)", "stroke-width": "4"
    }));

    container.appendChild(svg);
  }

  function renderAllCompositions() {
    Array.from(document.querySelectorAll(".albers-composition")).forEach(function (el, idx) {
      renderComposition(el, idx);
    });
  }

  /* ---------------------------------------------------------------------
   * Vignette title block: gather title / subtitle / author / date so they
   * lay out as one unit above the two-column body.
   * ------------------------------------------------------------------- */
  function initTitleBlock() {
    var body = document.body;
    var title = body.querySelector(":scope > h1.title");
    if (!title || body.querySelector(".albers-titleblock")) return;
    var header = document.createElement("header");
    header.className = "albers-titleblock";
    title.parentNode.insertBefore(header, title);
    header.appendChild(title);
    var sub = body.querySelector(":scope > h3.subtitle, :scope > p.subtitle");
    if (sub && sub.tagName === "H3") {
      // keep the outline h1 -> h2: the subtitle is prose, not a heading
      var p = document.createElement("p");
      p.className = sub.className;
      p.innerHTML = sub.innerHTML;
      sub.parentNode.replaceChild(p, sub);
      sub = p;
    }
    if (sub) header.appendChild(sub);
    var meta = Array.from(body.querySelectorAll(":scope > h4.author, :scope > h4.date, :scope > p.author, :scope > p.date"));
    if (meta.length) {
      var byline = document.createElement("div");
      byline.className = "albers-byline";
      meta.forEach(function (m) {
        // Demote to plain paragraphs: these are metadata, not headings.
        var p = document.createElement("p");
        p.className = m.className;
        p.innerHTML = m.innerHTML;
        byline.appendChild(p);
        m.parentNode.removeChild(m);
      });
      header.appendChild(byline);
    }
    if (!document.body.classList.contains("albers-no-plate")) {
      var plate = document.createElement("div");
      plate.className = "albers-plate-mark";
      plate.setAttribute("aria-hidden", "true");
      header.appendChild(plate);
    }
  }

  /* Dark twins: albers_vignette() renders each ggplot a second time with the
     night theme and emits it (hidden) right after the figure. Pair each twin
     with its figure, then show whichever matches the page theme. */
  /* Figures open at their natural size in a dialog (legible on phones). */
  function initFigureZoom() {
    if (typeof HTMLDialogElement === "undefined") return;
    var dlg = null;
    function open(src, alt, caption) {
      if (!dlg) {
        dlg = document.createElement("dialog");
        dlg.className = "albers-zoom";
        dlg.setAttribute("aria-label", "Figure, enlarged");
        dlg.addEventListener("click", function (e) { if (e.target === dlg) dlg.close(); });
        document.body.appendChild(dlg);
      }
      dlg.innerHTML = "";
      var size = document.createElement("button");
      size.type = "button";
      size.className = "albers-zoom-size";
      size.setAttribute("aria-pressed", "false");
      size.textContent = "Actual size";
      var close = document.createElement("button");
      close.type = "button";
      close.className = "albers-zoom-close";
      close.setAttribute("aria-label", "Close enlarged figure");
      close.textContent = "\u00d7";
      close.addEventListener("click", function (e) { e.stopPropagation(); dlg.close(); });
      dlg.appendChild(close);
      dlg.appendChild(size);
      var img = document.createElement("img");
      img.src = src;
      img.alt = alt || "";
      // at: the point (fractions of the image, and its place on screen) to
      // keep under the pointer; without one, start from the plot's centre
      // the caption band (at most 30vh) is fixed at the foot: reserve its
      // height below the plot, and again when the window changes size
      function reserveCaption() {
        var capEl = dlg.querySelector(".albers-zoom-stage > p");
        var on = img.classList.contains("is-actual");
        stage.style.paddingBottom = on && capEl ? (capEl.getBoundingClientRect().height + 16) + "px" : "";
      }
      if (!dlg.albersResize) {
        dlg.albersResize = true;
        window.addEventListener("resize", function () { if (dlg.open && dlg.albersReserve) dlg.albersReserve(); });
      }
      dlg.albersReserve = reserveCaption;
      function setActual(on, at) {
        img.classList.toggle("is-actual", on);
        dlg.classList.toggle("is-actual", on);
        size.setAttribute("aria-pressed", String(on));
        reserveCaption();
        if (!on) return;
        if (at) {
          var r = img.getBoundingClientRect();
          dlg.scrollLeft = r.left + dlg.scrollLeft + at.fx * r.width - at.cx;
          dlg.scrollTop = r.top + dlg.scrollTop + at.fy * r.height - at.cy;
        } else {
          dlg.scrollLeft = Math.max(0, (dlg.scrollWidth - dlg.clientWidth) / 2);
        }
      }
      size.addEventListener("click", function (e) { e.stopPropagation(); setActual(!img.classList.contains("is-actual")); });
      img.addEventListener("click", function (e) {
        if (img.classList.contains("is-actual")) { setActual(false); return; }
        var r = img.getBoundingClientRect();
        setActual(true, { fx: (e.clientX - r.left) / r.width, fy: (e.clientY - r.top) / r.height, cx: e.clientX, cy: e.clientY });
      });
      dlg.onkeydown = function (e) {
        if (e.key === "+" || e.key === "=") setActual(true);
        else if (e.key === "-" || e.key === "0") setActual(false);
        else if (img.classList.contains("is-actual")) {
          // the arrow keys pan the enlarged figure (Shift: further); Page
          // Up/Down a screen; Home/End to the top/bottom
          var step = e.shiftKey ? 200 : 40, page = dlg.clientHeight * 0.9;
          var d = { ArrowRight: [step, 0], ArrowLeft: [-step, 0], ArrowDown: [0, step], ArrowUp: [0, -step],
                    PageDown: [0, page], PageUp: [0, -page] }[e.key];
          if (d) { dlg.scrollBy(d[0], d[1]); e.preventDefault(); }
          else if (e.key === "Home") { dlg.scrollTo(0, 0); e.preventDefault(); }
          else if (e.key === "End") { dlg.scrollTo(dlg.scrollLeft, dlg.scrollHeight); e.preventDefault(); }
        }
      };
      var stage = document.createElement("div");
      stage.className = "albers-zoom-stage";
      stage.appendChild(img);
      dlg.appendChild(stage);
      var cap = document.createElement("p");
      // keep the caption's markup (its "Figure 1" label, emphasis, code)
      if (caption && caption.nodeType) {
        Array.from(caption.cloneNode(true).childNodes).forEach(function (n) { cap.appendChild(n); });
      }
      cap.id = "albers-zoom-caption";
      if (caption) { stage.appendChild(cap); dlg.setAttribute("aria-describedby", cap.id); }
      else dlg.removeAttribute("aria-describedby");
      dlg.showModal();
      close.focus();
    }
    Array.from(document.querySelectorAll("div.figure img, figure img, img.r-plt, p > img[src^='data:image'], img.albers-has-twin, img.albers-dark-twin")).forEach(function (img) {
      if (img.classList.contains("albers-zoomable")) return;
      if (img.closest("a")) return;
      img.classList.add("albers-zoomable");
      img.setAttribute("tabindex", "0");
      img.setAttribute("role", "button");
      img.setAttribute("aria-label", "Enlarge figure" + (img.alt ? ": " + img.alt : ""));
      function go() {
        var shown = img.closest("div.figure, figure") ?
          Array.from(img.closest("div.figure, figure").querySelectorAll("img")).filter(function (i) { return getComputedStyle(i).display !== "none"; })[0] || img : img;
        var capEl = img.closest("div.figure, figure");
        capEl = capEl && capEl.querySelector("p.caption, figcaption");
        open(shown.albersFull || shown.getAttribute("src"), shown.getAttribute("alt"), capEl && capEl.textContent.trim() ? capEl : null);
      }
      img.addEventListener("click", go);
      img.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); go(); } });
    });
  }

  // Code wraps like typeset code, not like prose. Each source line is cut into
  // segments that never break inside (white-space: pre); a line may break only
  // between segments, i.e. after ", ", after an opening "(", or after a spaced
  // binary operator -- never inside a name, a number or a string (strings and
  // comments are highlighted tokens, and only plain text between tokens is
  // cut). The text itself is untouched, so copy and find stay exact. A
  // segment too long for a phone may still wrap inside, as a last resort.
  // (a comparison is not a break point: "condition ==" / "\"incongruent\""
  // split the test from its value)
  var SEG_BREAK = /, +|\((?!\))|[)\]] +(?![-+*\/|%<>&=!~{])| (?:\+|-|\*|\/|\|>|%[^% ]*%|<-|->|&&|\|\||~|&|\|) +/g;
  var SEG_LONG = 28;

  function segmentLine(line) {
    var kids = Array.from(line.childNodes);
    var start = 0;
    // pandoc's (empty) line anchor stays first, outside the segments
    while (start < kids.length && kids[start].nodeType === 1 && kids[start].tagName === "A" && !kids[start].textContent) start++;
    var nodes = kids.slice(start);
    // the indentation holds on to the first character (CSS adds a word
    // joiner after it), so a line never wraps right after its indent
    // (YAML's highlighter puts the indent in its own span, or at the start
    // of the first token's: then that token is the no-wrap lead)
    if (nodes.length && nodes[0].nodeType === 1 && /^ +$/.test(nodes[0].textContent)) {
      nodes[0].classList.add("albers-indent");
    } else if (nodes.length && nodes[0].nodeType === 1 && /^ +\S/.test(nodes[0].textContent) && nodes[0].textContent.length <= 36) {
      nodes[0].classList.add("albers-lead");
      // its leading spaces are indentation too (drawn at half width on phones)
      var t0 = nodes[0].firstChild;
      if (t0 && t0.nodeType === 3) {
        var sp = /^ +/.exec(t0.nodeValue);
        if (sp && t0.nodeValue.length > sp[0].length) {
          t0.splitText(sp[0].length);
          var ind0 = document.createElement("span");
          ind0.className = "albers-indent";
          nodes[0].insertBefore(ind0, t0);
          ind0.appendChild(t0);
        }
      }
    } else if (nodes.length && nodes[0].nodeType === 3) {
      var ind = /^ +/.exec(nodes[0].nodeValue);
      if (ind) {
        if (nodes[0].nodeValue.length > ind[0].length) nodes[0].splitText(ind[0].length);
        var indSpan = document.createElement("span");
        indSpan.className = "albers-indent";
        line.insertBefore(indSpan, nodes[0]);
        indSpan.appendChild(nodes[0]);
        kids = Array.from(line.childNodes);
        nodes = kids.slice(start);
      }
    }
    // the indent and the first token never part (one no-wrap lead)
    var first = nodes[0];
    if (first && first.nodeType === 1 && first.classList.contains("albers-indent") && first.nextSibling) {
      var tok = first.nextSibling;
      if (tok.nodeType === 3) {
        var w = /^\S+/.exec(tok.nodeValue);
        if (w && tok.nodeValue.length > w[0].length) tok.splitText(w[0].length);
      }
      if (tok.textContent.length <= 36) {
        var lead = document.createElement("span");
        lead.className = "albers-lead";
        line.insertBefore(lead, first);
        lead.appendChild(first);
        lead.appendChild(tok);
        kids = Array.from(line.childNodes);
        nodes = kids.slice(start);
      }
    }
    var text = nodes.map(function (n) { return n.textContent; }).join("");
    if (text.length <= SEG_LONG) {
      // a short line is one unbreakable segment
      if (nodes.length && /\S/.test(text)) wrapSeg(nodes, text.length);
      return;
    }
    var lead = /^\s*/.exec(text)[0].length;
    // the continuation hangs just after the line's first open bracket (so a
    // wrapped argument lines up under the first one, as it would be written)
    // -- only when that bracket is still open at the line's end, or closes
    // just before it: after a call that closed mid-line ("ggplot(d) + ..."),
    // a continuation is not one of its arguments and takes the plain hang
    var bracket = -1, close = -1, depth = 0, quote = null;
    for (var qi = lead; qi < text.length; qi++) {
      var qc = text.charAt(qi);
      // (a backslash escapes the next character: "\\" closes its string)
      if (quote) { if (qc === "\\") qi++; else if (qc === quote) quote = null; continue; }
      if (qc === '"' || qc === "'" || qc === "`") { quote = qc; continue; }
      if (qc === "#") break;
      if (qc === "(" || qc === "[") { if (bracket < 0) bracket = qi; depth++; }
      else if ((qc === ")" || qc === "]") && bracket >= 0 && --depth === 0) { close = qi; break; }
    }
    var tail = close < 0 ? "" : text.slice(close + 1).replace(/\s*#.*$/, "").trim();
    if (bracket >= 0 && bracket + 1 > lead && /^[)\],;]*\s*(?:\{|\+|\|>|%>%|,)?$/.test(tail)) {
      line.setAttribute("data-hang", String(bracket + 1));
    }
    var offsets = [], m;
    SEG_BREAK.lastIndex = 0;
    // a short bracketed call (it closes within 24 characters) is kept whole:
    // no break after its "(" or inside it, so breaks fall at shallower points
    // (no "seq_len(" / "n_trial)," or "round(" / "agg$rt, 1))" stubs)
    var shut = [], stack = [];
    for (var ci = 0; ci < text.length; ci++) {
      var ch = text.charAt(ci);
      if (ch === "(" || ch === "[") stack.push(ci);
      else if ((ch === ")" || ch === "]") && stack.length) {
        var open = stack.pop();
        if (ci - open <= 24) shut.push([open, ci]);
      }
    }
    function inShort(o) {
      for (var k = 0; k < shut.length; k++) if (o > shut[k][0] && o <= shut[k][1]) return true;
      return false;
    }
    // string literals (with escapes) are never split at a segment break
    var strs = [], sq = null, s0 = 0;
    for (var si = 0; si < text.length; si++) {
      var sc = text.charAt(si);
      if (sq) { if (sc === "\\") si++; else if (sc === sq) { strs.push([s0, si + 1]); sq = null; } }
      else if (sc === '"' || sc === "'") { sq = sc; s0 = si; }
      else if (sc === "#") break;
    }
    function inString(o) {
      for (var k = 0; k < strs.length; k++) if (o > strs[k][0] && o < strs[k][1]) return true;
      return false;
    }
    var lambdas = [], lm, LAM = /\\\([^()]*\)/g;
    // (a "\\(" inside a string literal is no lambda)
    while ((lm = LAM.exec(text))) if (!inString(lm.index + 1)) lambdas.push([lm.index, lm.index + lm[0].length]);
    function inLambda(o) {
      for (var k = 0; k < lambdas.length; k++) if (o > lambdas[k][0] && o < lambdas[k][1]) return true;
      return false;
    }
    var soft = [];
    while ((m = SEG_BREAK.exec(text))) {
      var o = m.index + m[0].length;
      if (!(o > lead && o < text.length)) continue;
      // never inside R's lambda head "\(x, y)" (it stays one text node,
      // for glueLambda)
      if (text.charAt(o - 1) === "\\" || inLambda(o) || inString(o)) continue;
      // inside a short call: not a segment boundary, but a soft break that
      // counts only if its segment must wrap (and beats one after "::")
      if (inShort(o)) soft.push(o); else offsets.push(o);
    }
    if (!offsets.length && !soft.length) {
      // no legal break at all (a long comment, name or string): one segment
      // that may wrap inside, as the last resort
      wrapSeg(nodes, text.length);
      return;
    }
    // split plain text nodes at the offsets (soft ones too); an offset inside
    // a token is dropped
    var all = offsets.concat(soft).sort(function (a, b) { return a - b; });
    var isSoft = {};
    soft.forEach(function (o) { isSoft[o] = true; });
    var parts = [], pos = 0, oi = 0;
    nodes.forEach(function (n) {
      var len = n.textContent.length;
      if (n.nodeType === 3) {
        while (oi < all.length && all[oi] <= pos) oi++;
        while (oi < all.length && all[oi] < pos + len) {
          var local = all[oi] - pos;
          if (local > 0) {
            var rest = n.splitText(local);
            parts.push(n);
            pos += local;
            len -= local;
            n = rest;
          }
          oi++;
        }
      }
      parts.push(n);
      pos += len;
    });
    // group into segments at the offsets
    var cut = {};
    offsets.forEach(function (o) { cut[o] = true; });
    var seg = [], at = 0, segs = [];
    parts.forEach(function (n) {
      seg.push(n);
      at += n.textContent.length;
      if (cut[at]) { segs.push(seg); seg = []; }
      else if (isSoft[at] && n.parentNode) {
        var w = document.createElement("wbr");
        n.parentNode.insertBefore(w, n.nextSibling);
        seg.push(w);
      }
    });
    if (seg.length) segs.push(seg);
    if (segs.length < 2) { wrapSeg([].concat.apply([], segs), text.length); return; }
    segs.forEach(function (g, i) {
      var span = wrapSeg(g, g.map(function (n) { return n.textContent; }).join("").length);
      if (i < segs.length - 1) span.parentNode.insertBefore(document.createElement("wbr"), span.nextSibling);
    });
  }

  function wrapSeg(nodes, len) {
    var span = document.createElement("span");
    span.className = len > SEG_LONG ? "albers-seg albers-seg-long" : "albers-seg";
    nodes[0].parentNode.insertBefore(span, nodes[0]);
    nodes.forEach(function (n) { span.appendChild(n); });
    // any segment can have to wrap on a narrow screen: its units are glued
    stringRuns(span); glueOperators(span); glueBrace(span); glueAtValue(span); glueLambda(span); softBreaks(span);
    if (len > SEG_LONG) pathBreaks(span);
    // an argument name (.at) stays whole only when it is short: a long YAML
    // value or name may wrap rather than run off the block
    Array.from(span.querySelectorAll(".at")).forEach(function (at) {
      if (at.textContent.length <= 24) at.classList.add("albers-at-short");
    });
    glueNames(span);
    return span;
  }

  // a string and the punctuation that closes it ('"...")', '"...", ') as
  // one inline run: when its segment has to wrap, the run is kept whole if it
  // fits the line's room (fitSegments), so the break falls before the string
  // does this text hold a whole string literal (its closing quote, not an
  // escaped one)?
  function literalClosed(t) {
    var q = t.charAt(0);
    if (q !== '"' && q !== "'") return true;
    for (var i = 1; i < t.length; i++) {
      if (t.charAt(i) === "\\") i++;
      else if (t.charAt(i) === q) return true;
    }
    return false;
  }
  function stringRuns(span) {
    // (pandoc writes an escape inside a string as its own .sc span between
    // .st spans -- "\\(" is .st .sc .st -- so a literal's adjoining pieces
    // are joined into one first; the line breaker allowed "\\" / "(")
    Array.from(span.querySelectorAll(".st")).forEach(function (st) {
      if (!st.parentNode) return;
      var nx = st.nextSibling;
      while (nx && nx.nodeType === 1 && (nx.classList.contains("sc") || nx.classList.contains("st")) && !literalClosed(st.textContent)) {
        var after = nx.nextSibling;
        st.appendChild(nx);
        nx = after;
      }
    });
    Array.from(span.querySelectorAll(".st")).forEach(function (st) {
      if (st.parentNode && st.parentNode.closest && st.parentNode.closest(".st")) return;
      var run = [st], next = st.nextSibling;
      // downlit puts ")" in its own .op span
      while (next && next.nodeType === 1 && /^[)\]}]+$/.test(next.textContent)) { run.push(next); next = next.nextSibling; }
      if (next && next.nodeType === 3) {
        var m = /^[)\]}]*(?:, ?)?/.exec(next.nodeValue);
        if (m && m[0]) {
          if (next.nodeValue.length > m[0].length) next.splitText(m[0].length);
          run.push(next);
        }
      }
      // an operator right after the string joins it: '"...") +' never
      // leaves "+" alone on a line
      var after = run[run.length - 1].nextSibling;
      if (after && after.nodeType === 3 && /^ +$/.test(after.nodeValue) && after.nextSibling &&
          after.nextSibling.nodeType === 1 && OPS.test(after.nextSibling.textContent)) {
        run.push(after, after.nextSibling);
      }
      var u = document.createElement("span");
      u.className = "albers-strq" + (st.textContent.length > 30 ? " albers-strq-long" : "");
      st.parentNode.insertBefore(u, st);
      u.appendChild(st);
      // what closes the string ('") +') never breaks, even when a long
      // string itself has to
      if (run.length > 1) {
        var tail = document.createElement("span");
        tail.className = "albers-tail";
        u.appendChild(tail);
        run.slice(1).forEach(function (n) { tail.appendChild(n); });
      }
    });
  }

  // a spaced binary operator stays with the token before it, so no line
  // starts with "+", "<-" or "|>"
  // characters wider (or narrower) than one monospace column
  var WIDE = /[^\u0000-\u024f\u2000-\u206f\u2190-\u22ff]/;
  var OPS = /^(\+|-|\*|\/|\|>|%[^% ]*%|<-|->|==|!=|<=|>=|<|>|&&|\|\||~|&|\|)$/;
  // a comparison also keeps its value ("rt > 200", never "rt >" / "200")
  var CMP = /^(==|!=|<=|>=|<|>)$/;
  function glueOperators(span) {
    Array.from(span.querySelectorAll("span")).forEach(function (op) {
      if (!OPS.test(op.textContent) || op.parentNode !== span) return;
      var ws = op.previousSibling, tok;
      if (!ws || ws.nodeType !== 3 || !/ +$/.test(ws.nodeValue)) return;
      if (!/^ +$/.test(ws.nodeValue)) {
        // "x) " before the operator: its trailing spaces are the gap
        ws = ws.splitText(ws.nodeValue.search(/ +$/));
      }
      tok = ws.previousSibling;
      if (!tok) return;
      if (tok.nodeType === 3) {
        var m = /\S+$/.exec(tok.nodeValue);
        if (!m) return;
        if (m.index > 0) tok = tok.splitText(m.index);
      }
      if (tok.textContent.length > 24) {
        // a long operand stays apart, but the operator keeps the operand's
        // closing punctuation (or its last few characters), so it never
        // starts a row alone ("xxx))" / "+"); a comparison keeps its value
        if (tok.nodeType === 3) {
          var tail = /(?:[)\]}]+|\w{1,8})$/.exec(tok.nodeValue);
          if (tail) {
            var piece = tok.splitText(tail.index);
            var t2 = document.createElement("span");
            t2.className = "albers-opq";
            span.insertBefore(t2, piece);
            t2.appendChild(piece); t2.appendChild(ws); t2.appendChild(op);
            if (CMP.test(op.textContent)) glueValue(span, t2);
            return;
          }
        }
        if (CMP.test(op.textContent)) glueValue(span, op);
        return;
      }
      var g = document.createElement("span");
      g.className = "albers-opq";
      span.insertBefore(g, tok);
      g.appendChild(tok); g.appendChild(ws); g.appendChild(op);
      if (CMP.test(op.textContent)) glueValue(span, g);
    });
  }

  // pandoc's `name =` (.at) stays with a short value after it, so the
  // wrap never leaves "trial =" at a line end with "seq_len(" below
  function glueAtValue(span) {
    Array.from(span.querySelectorAll(":scope > .at")).forEach(function (at) {
      if (!/=\s*$/.test(at.textContent)) return;
      glueValue(span, at);
    });
  }

  // glue `head`, the space after it and the first word of what follows (a
  // short value) into one no-break unit
  function glueValue(span, head) {
    var ws = head.nextSibling, val;
    if (!ws || ws.nodeType !== 3 || !/^ +/.test(ws.nodeValue)) return;
    if (/^ +$/.test(ws.nodeValue)) {
      val = ws.nextSibling;
      if (!val) return;
      if (val.nodeType === 3) {
        var w = /^\S+/.exec(val.nodeValue);
        if (!w) return;
        if (val.nodeValue.length > w[0].length) val.splitText(w[0].length);
      }
    } else {
      // " agg..." : the space and the word share one text node
      var m = /^( +)(\S+)/.exec(ws.nodeValue);
      if (ws.nodeValue.length > m[0].length) ws.splitText(m[0].length);
      val = null;
    }
    if ((val ? val.textContent.length : 0) > 16) return;
    var g = document.createElement("span");
    g.className = "albers-nameq";
    span.insertBefore(g, head);
    g.appendChild(head); g.appendChild(ws); if (val) g.appendChild(val);
  }

  // R's lambda head "\(x, y)" never breaks, on screen or paper (the
  // browser's own line breaker allowed "\" / "(x)")
  function glueLambda(span) {
    var tw = document.createTreeWalker(span, NodeFilter.SHOW_TEXT), t, hits = [];
    while ((t = tw.nextNode())) if (/\\\(/.test(t.nodeValue)) hits.push(t);
    hits.forEach(function (node) {
      if (node.parentNode.closest && node.parentNode.closest(".st, .co")) return;
      var m = /\\\([^()]*\)/.exec(node.nodeValue);
      if (!m) return;
      var head = m.index > 0 ? node.splitText(m.index) : node;
      if (head.nodeValue.length > m[0].length) head.splitText(m[0].length);
      var g = document.createElement("span");
      g.className = "albers-lam";
      head.parentNode.insertBefore(g, head);
      g.appendChild(head);
    });
  }

  // an opening brace at the end of a line stays with what comes before it
  // ("function(x) {", "if (x) {"), never alone on a line
  function glueBrace(span) {
    var last = span.lastChild;
    while (last && last.nodeType === 3 && !/\S/.test(last.nodeValue)) last = last.previousSibling;
    if (!last || last.nodeType !== 3) return;
    var m = /(\S*\s+\{\s*)$/.exec(last.nodeValue);
    if (!m) return;
    var piece = m.index > 0 ? last.splitText(m.index) : last;
    var run = [piece];
    // "( ... ) {": take the token before too, if the text started with the space
    if (/^\s/.test(piece.nodeValue) && piece.previousSibling) run.unshift(piece.previousSibling);
    var g = document.createElement("span");
    g.className = "albers-opq";
    span.insertBefore(g, run[0]);
    run.forEach(function (n) { g.appendChild(n); });
  }

  // a long string (a path, a URL) may break after "/", "_", "-" or "."
  function pathBreaks(span) {
    Array.from(span.querySelectorAll(".albers-strq-long .st")).forEach(function (st) {
      var t = st.firstChild;
      if (!t || t.nodeType !== 3 || st.childNodes.length !== 1) return;
      var parts = t.nodeValue.replace(/([\/_.\-?&=])/g, "$1\u0000").split("\u0000");
      if (parts.length < 2) return;
      st.textContent = "";
      parts.forEach(function (p, i) {
        if (!p) return;
        if (i) st.appendChild(document.createElement("wbr"));
        st.appendChild(document.createTextNode(p));
      });
    });
  }

  // downlit writes `name <span class="op">=</span>`: the name and its "="
  // never part, so "=" cannot start a line even inside a wrapping segment
  function glueNames(span) {
    Array.from(span.querySelectorAll(".op")).forEach(function (eq) {
      if (eq.textContent !== "=") return;
      var prev = eq.previousSibling;
      if (!prev || prev.nodeType !== 3) return;
      var m = /(\S+\s*)$/.exec(prev.nodeValue);
      if (!m) return;
      var name = m.index > 0 ? prev.splitText(m.index) : prev;
      var g = document.createElement("span");
      g.className = "albers-nameq";
      eq.parentNode.insertBefore(g, name);
      g.appendChild(name);
      g.appendChild(eq);
      if (g.parentNode === span) glueValue(span, g);
    });
  }

  // inside a long segment, prefer breaking after "::" and "$" (pkg::fun,
  // df$col) to breaking inside a name
  // (they are ranked below every other break: see nsNeeded)
  function softBreaks(span) {
    var tw = document.createTreeWalker(span, NodeFilter.SHOW_TEXT);
    var texts = [], t;
    while ((t = tw.nextNode())) texts.push(t);
    texts.forEach(function (node) {
      // never inside a string or comment token
      if (node.parentNode && node.parentNode.closest && node.parentNode.closest(".st, .co, .albers-strq")) return;
      var m, re = /::|\$/g, cuts = [];
      while ((m = re.exec(node.nodeValue))) {
        var end = m.index + m[0].length;
        if (end < node.nodeValue.length || node.parentNode !== span) cuts.push(end);
      }
      for (var i = cuts.length - 1; i >= 0; i--) {
        var after = cuts[i] < node.nodeValue.length ? node.splitText(cuts[i]) : node.nextSibling;
        var wbr = document.createElement("wbr");
        // the last resort: shown only when no other break can make the
        // segment fit (fitSegments), so "dplyr::" never ends a line that
        // could have broken after "(" or "|>"
        wbr.className = "albers-wbr-ns";
        node.parentNode.insertBefore(wbr, after);
      }
    });
  }

  // Only a segment wider than the room its line has (the column less the
  // hang, where a continuation starts) may wrap inside. Measured after layout
  // and again on resize, when fonts arrive, and before printing.
  // a segment's width without its trailing spaces (they hang at a line end
  // and need no room)
  function inkWidth(el) {
    var tw = document.createTreeWalker(el, NodeFilter.SHOW_TEXT), t, texts = [];
    while ((t = tw.nextNode())) texts.push(t);
    var last = -1;
    texts.forEach(function (n, i) { if (/\S/.test(n.nodeValue)) last = i; });
    if (last < 0) return el.getBoundingClientRect().width;
    // the sum of the glyph boxes of its text (not of its elements), so a
    // segment that is wrapping now still measures as one line
    var r = document.createRange(), sum = 0;
    texts.slice(0, last + 1).forEach(function (n, i) {
      r.setStart(n, 0);
      r.setEnd(n, i === last ? n.nodeValue.replace(/\s+$/, "").length : n.nodeValue.length);
      Array.from(r.getClientRects()).forEach(function (q) { sum += q.width; });
    });
    return sum;
  }


  // one monospace character's width in a code element (measured once per
  // font, per refit)
  var charWidths = {};
  function charWidth(code) {
    var cs = getComputedStyle(code), key = cs.fontFamily + "|" + cs.fontSize + "|" + cs.fontWeight + "|" + cs.letterSpacing;
    if (charWidths[key]) return charWidths[key];
    var w = measureChar(code);
    if (w) charWidths[key] = w;
    return w;
  }
  function measureChar(code) {
    // a probe in the block's font, outside any block (a block's first
    // character may be hidden or not laid out yet)
    // (it stays in the page, so measuring again changes nothing to lay out)
    var cs = getComputedStyle(code);
    if (!charProbe) {
      // inside a fixed, empty, clipped box: it adds no scroll in any
      // direction (a left-offset probe gave right-to-left pages 10,000 px)
      var box = document.createElement("div");
      box.setAttribute("aria-hidden", "true");
      box.className = "albers-char-probe";
      box.style.cssText = "position:fixed;top:0;left:0;width:0;height:0;overflow:hidden;visibility:hidden;pointer-events:none;contain:strict";
      charProbe = document.createElement("span");
      charProbe.textContent = "0000000000";
      box.appendChild(charProbe);
      document.body.appendChild(box);
    }
    var css = "position:absolute;white-space:pre;top:0;left:0;font-family:" + cs.fontFamily +
      ";font-size:" + cs.fontSize + ";font-weight:" + cs.fontWeight + ";letter-spacing:" + cs.letterSpacing +
      ";font-feature-settings:" + cs.fontFeatureSettings + ";font-variant-ligatures:" + cs.fontVariantLigatures;
    // (compared as written: the browser normalizes cssText)
    if (charProbe.getAttribute("data-font") !== css) { charProbe.style.cssText = css; charProbe.setAttribute("data-font", css); }
    return charProbe.getBoundingClientRect().width / 10;
  }
  var charProbe = null, printingCode = false;


  // (`lines`, with a single code element: fit only those lines, the ones a
  // chunk of preparation just added)
  // the width a block's lines have (fractional: clientWidth rounds, and a
  // segment that exactly fills its line must not be judged too wide)
  function roomOf(pre) {
    var ps = getComputedStyle(pre);
    return pre.getBoundingClientRect().width - parseFloat(ps.borderLeftWidth) - parseFloat(ps.borderRightWidth) -
      parseFloat(ps.paddingLeft) - parseFloat(ps.paddingRight) - (pre.offsetWidth - pre.clientWidth - parseFloat(ps.borderLeftWidth) - parseFloat(ps.borderRightWidth));
  }

  function fitSegments(only, lines) {
    // batched: all reads (one layout), then only the class changes -- no
    // clearing pass first: it forced a second layout of the whole page on
    // every resize (seconds on a long listing). Widths are arithmetic (or
    // summed glyph boxes), so a segment's current wrapping does not matter.
    var codes = only ? [].concat(only) : Array.from(document.querySelectorAll("pre > code.sourceCode[data-albers-seg]"));
    var jobs = [], stale = [];
    codes.forEach(function (code) {
      (lines || [code]).forEach(function (root) {
        Array.from(root.querySelectorAll(".is-wrap, .is-keep, .is-on")).forEach(function (el) { stale.push(el); });
      });
    });
    codes.forEach(function (code) {
      var room = roomOf(code.parentElement);
      if (!(room > 0)) return;
      // code is monospace: a line whose characters cannot fill the room
      // needs no measuring (most lines), so only long ones are read
      var cw = charWidth(code) || 8;
      // widths by arithmetic: characters x one character's width
      // (indentation included: it is drawn at full width everywhere)
      var iw = cw;
      // (paper takes the author's own indentation, whatever window printed it)
      var phone = !printingCode && window.matchMedia && window.matchMedia("(max-width: 620px)").matches;
      function width(el) {
        var text = el.textContent.replace(/\s+$/, "");
        // wide scripts (CJK, emoji) are not one character wide: measure them
        if (WIDE.test(text)) return inkWidth(el);
        var ind = 0;
        Array.from(el.querySelectorAll(".albers-indent")).forEach(function (s) { ind += s.textContent.length; });
        // (indentation at the width it is drawn: on a phone an author's
        // alignment may be redrawn at the call's depth, below)
        var ln = ind && el.closest ? el.closest("pre > code > span") : null;
        var shown = ln && drawn.has(ln) ? drawn.get(ln) : ind;
        return (text.length - ind) * cw + shown * iw;
      }
      // On a phone an author's alignment under an open bracket is redrawn at
      // the call's depth: the aligned lines, and everything nested in them,
      // move left to the opener's indent + 2ch -- the column its own wrapped
      // arguments take -- so siblings share one column and children stay
      // deeper than their parent. Only the drawn width of the indentation
      // changes (letter-spacing), never the text. Arithmetic only.
      var drawn = new Map(), redrawnOpeners = new Set();
      var all = Array.from(code.querySelectorAll(":scope > span:not(.albers-out-block)"));
      if (phone) {
        var orig = all.map(function (ln) { return parseInt(ln.getAttribute("data-ind") || "0", 10) || 0; });
        var cur = orig.slice();
        var blank = all.map(function (ln) { return !/\S/.test(ln.textContent); });
        all.forEach(function (ln, i) {
          var at = parseInt(ln.getAttribute("data-hang") || "", 10);
          if (!(at > 0)) return;
          var j = i + 1;
          while (j < all.length && blank[j]) j++;
          if (j >= all.length || orig[j] !== at) return;
          // the bracket's column as drawn, and the column it moves to
          var extra = (at - (orig[i] - cur[i])) - (cur[i] + 2);
          if (extra <= 0) return;
          redrawnOpeners.add(ln);
          for (var k = j; k < all.length && (blank[k] || orig[k] >= at); k++) if (!blank[k]) cur[k] -= extra;
        });
        all.forEach(function (ln, i) { if (cur[i] !== orig[i]) drawn.set(ln, cur[i]); });
      }
      all.forEach(function (ln) {
        var o = parseInt(ln.getAttribute("data-ind") || "0", 10) || 0;
        if (drawn.has(ln)) {
          var d = drawn.get(ln);
          // (over the indent's o spaces: its zero-width joiner takes no
          // letter-spacing)
          jobs.push({ el: ln, prop: "--ind", val: d + "ch" });
          jobs.push({ el: ln, prop: "--ind-ls", val: ((d - o) / o).toFixed(4) + "ch" });
        } else if (ln.style.getPropertyValue("--ind-ls")) {
          jobs.push({ el: ln, prop: "--ind", val: o ? o + "ch" : null });
          jobs.push({ el: ln, prop: "--ind-ls", val: null });
        }
      });
      // the "::" and "$" breaks inside a run (between the segment's ordinary
      // break points: spaces outside glued units, visible <wbr>s) that is
      // wider than the room -- only those may be used
      function nsNeeded(seg, fit) {
        var run = 0, inRun = [], on = [];
        function end() { if (run * cw > fit) on.push.apply(on, inRun); run = 0; inRun = []; }
        (function walk(node, glued) {
          Array.from(node.childNodes).forEach(function (c) {
            if (c.nodeType === 3) {
              if (glued) { run += c.nodeValue.length; return; }
              c.nodeValue.split(/( +)/).forEach(function (p) {
                if (/^ +$/.test(p)) end(); else run += p.length;
              });
            } else if (c.tagName === "WBR") {
              if (c.classList.contains("albers-wbr-ns")) inRun.push(c);
              else if (!glued) end();
            } else {
              walk(c, glued || c.matches(".albers-nameq, .albers-opq, .albers-lead, .albers-indent, .albers-lam, .albers-tail, .albers-strq:not(.albers-strq-long)"));
            }
          });
        })(seg, false);
        end();
        return on;
      }
      (lines || Array.from(code.querySelectorAll(":scope > span:not(.albers-out-block)"))).forEach(function (line) {
        var lt = line.textContent.replace(/\s+$/, "");
        // (wide scripts are wider than their count: those lines are measured)
        if (lt.length * cw <= room - 2 * cw && !WIDE.test(lt)) return;
        var segs = line.querySelectorAll(":scope > .albers-seg");
        if (!segs.length) return;
        var key = String(drawn.has(line) ? drawn.get(line) : (parseInt(line.getAttribute("data-ind") || "0", 10) || 0));
        // rows the line takes with continuations at hang h (greedy, as the
        // browser fills them: the first row has the whole room)
        var segW = null;
        // (a segment's trailing spaces take room when another follows on
        // its row: without them the fill was predicted too loose, and a hang
        // said to cost no rows added 1,000 on a long page)
        var segFull = null;
        function rowsAt(h) {
          if (!segW) {
            segW = Array.from(segs).map(width);
            segFull = Array.from(segs).map(function (sg, i) { return segW[i] + (sg.textContent.length - sg.textContent.replace(/\s+$/, "").length) * cw; });
          }
          var rows = 1, cap = room, cur = 0, next = Math.max(cw, room - h);
          segW.forEach(function (w, i) {
            if (cur + w <= cap + 0.5) { cur += segFull[i]; return; }
            if (cur > 0) { rows++; cap = next; cur = 0; }
            while (w > cap + 0.5) { rows++; w -= cap; cap = next; }
            cur = w + (segFull[i] - segW[i]);
          });
          return rows;
        }
        // the plain hang, as the CSS draws it: the indent plus 4ch, or plus
        // 2ch on phones (where indentation is drawn at full width too: halved,
        // siblings aligned by their author landed in different columns)
        var hang = (parseInt(key, 10) || 0) * iw + (phone ? 2 : 4) * cw;
        var plainHang = hang;
        var at = parseInt(line.getAttribute("data-hang") || "", 10);
        // (an opener whose aligned lines were redrawn at its depth wraps at
        // that same column: the plain hang)
        if (at > 0 && redrawnOpeners.has(line)) jobs.push({ el: line, hang: null });
        else if (at > 0) {
          // capped at half the room (a quarter on phones) so deep calls
          // still leave space to wrap into
          var indN = parseInt(key, 10) || 0;
          // (never left of the plain hang: a deeply indented line's
          // continuation went back under its own indent)
          // (an author who aligned the next line under this bracket set the
          // column: continuations take it too, so siblings share one)
          var nxt = line.nextElementSibling;
          while (nxt && nxt.classList.contains("albers-out-block")) nxt = nxt.nextElementSibling;
          var aligned = !!nxt && (parseInt(nxt.style.getPropertyValue("--ind"), 10) || 0) === at;
          var want = Math.max(0, at - indN) * cw + Math.min(indN, at) * iw;
          // (the author's column is honoured on a phone while every segment
          // fits beside it, or at least 20 columns remain there to wrap in;
          // on a narrower screen the usual rules decide)
          if (aligned && phone && room - want < 20 * cw) {
            rowsAt(0);
            if (segW.some(function (w) { return w > room - want + 0.5; })) aligned = false;
          }
          var cap0 = room * (phone && !aligned ? 0.25 : 0.5), cap = Math.max(cap0, hang);
          var px = Math.min(want, cap);
          // ... unless a string or pair then fits nowhere: on a very narrow
          // column the smaller hang beats breaking the unit
          var units = Array.from(line.querySelectorAll(".albers-strq, .albers-nameq")).map(width);
          var px0 = Math.min(want, cap0);
          if (px0 < px && units.some(function (w) { return w > room - px + 0.5 && w <= room - px0 + 0.5; })) px = px0;
          // but a string or a `name = value` pair that fits at the plain
          // hang and not at the bracket is worth more than the alignment:
          // keep the plain hang then
          var squeezed = Array.from(line.querySelectorAll(".albers-strq, .albers-nameq")).some(function (q) {
            var w = width(q);
            return w > room - px + 0.5 && w <= room - hang + 0.5;
          });
          // on whole character columns (10.3 ch aligned with nothing)
          px = Math.round(px / cw) * cw;
          // on a phone the alignment is taken only when it costs no rows
          // (one argument per row made code half as long again)
          // (whether the alignment costs rows is measured after the writes
          // below: an estimate of the browser's fill added 1,000 rows)
          var check = phone && !aligned && px > hang;
          if (px > 0 && !squeezed) { hang = px; jobs.push({ el: line, hang: px, check: check }); }
          else jobs.push({ el: line, hang: null });
        }
        Array.from(line.querySelectorAll(".albers-lead")).forEach(function (lead) {
          jobs.push({ el: lead, cls: "is-wrap", over: width(lead) > room + 0.5 });
        });
        // the classes for continuations at hang h (also computed for the
        // plain hang when a phone may take it back, below)
        function segJobs(h, out) {
          var fit = room - h + 0.5;
          Array.from(segs).forEach(function (seg) {
            var over = width(seg) > fit;
            out.push({ el: seg, cls: "is-wrap", over: over });
            if (!over) return;
            Array.from(seg.querySelectorAll(".albers-strq")).forEach(function (q) {
              out.push({ el: q, cls: "is-keep", over: width(q) <= fit });
            });
            // a glued unit ("name = value", "x +", "{") wider than the room
            // wraps inside rather than run under the edge
            Array.from(seg.querySelectorAll(".albers-nameq, .albers-opq")).forEach(function (g) {
              out.push({ el: g, cls: "is-wrap", over: width(g) > fit });
            });
            // the breaks after "::" and "$" only when no other can make it fit
            if (seg.querySelector(".albers-wbr-ns")) {
              nsNeeded(seg, fit).forEach(function (w) { out.push({ el: w, cls: "is-on", over: true }); });
            }
          });
        }
        var from = jobs.length;
        segJobs(hang, jobs);
        var hangJob = jobs.filter(function (j) { return j.el === line && j.check; })[0];
        if (hangJob) { hangJob.main = jobs.slice(from); hangJob.alt = []; segJobs(plainHang, hangJob.alt); }
      });
    });
    var want = new Map();
    jobs.forEach(function (j) {
      if (j.prop) {
        if (j.val === null) j.el.style.removeProperty(j.prop);
        else if (j.el.style.getPropertyValue(j.prop) !== j.val) j.el.style.setProperty(j.prop, j.val);
      }
      else if (j.hang === null) { if (j.el.style.getPropertyValue("--hang")) j.el.style.removeProperty("--hang"); }
      else if (j.hang !== undefined) { if (j.el.style.getPropertyValue("--hang") !== j.hang + "px") j.el.style.setProperty("--hang", j.hang + "px"); }
      else if (j.over) { var w = want.get(j.el) || []; w.push(j.cls); want.set(j.el, w); }
    });
    // classes no longer wanted come off; wanted ones go on (a class that is
    // already right is not touched, so an unchanged block is not re-laid out)
    stale.forEach(function (el) {
      var w = want.get(el) || [];
      ["is-wrap", "is-keep", "is-on"].forEach(function (c) { if (w.indexOf(c) < 0 && el.classList.contains(c)) el.classList.remove(c); });
    });
    want.forEach(function (cls, el) {
      cls.forEach(function (c) { if (!el.classList.contains(c)) el.classList.add(c); });
    });
    // On a phone a bracket alignment is kept only if the line is no taller
    // with it than with the plain hang: both are measured (two layouts, all
    // lines at once), not predicted.
    var cand = jobs.filter(function (j) { return j.check; });
    if (cand.length) {
      function wrapAs(j, list) {
        Array.from(j.el.querySelectorAll(".is-wrap, .is-keep, .is-on")).forEach(function (el) { el.classList.remove("is-wrap", "is-keep", "is-on"); });
        (list || []).forEach(function (a) { if (a.over) a.el.classList.add(a.cls); });
      }
      var withHang = cand.map(function (j) { return j.el.getBoundingClientRect().height; });
      // (each at the plain hang with the wrapping computed for it)
      cand.forEach(function (j) { j.el.style.removeProperty("--hang"); wrapAs(j, j.alt); });
      var plain = cand.map(function (j) { return j.el.getBoundingClientRect().height; });
      cand.forEach(function (j, i) {
        if (withHang[i] <= plain[i] + 0.5) { j.el.style.setProperty("--hang", j.hang + "px"); wrapAs(j, j.main); }
      });
    }
  }

  function watchSegments() {
    // one refit per frame, whatever asks for it (a callback never passes its
    // event on to fitSegments, which takes an optional code element)
    var queued = false;
    // a page opened at a #section keeps that section in place while code
    // above it is refitted (fonts arriving change the wrapping), until the
    // reader scrolls on their own
    // -- any scroll that moves it away (a screen reader, the find bar)
    // counts as the reader's
    // (a reload or Back/Forward of a long page lands on the place albers.js
    // saved; of a short one, the browser restores its own offset)
    var nav = performance.getEntriesByType ? performance.getEntriesByType("navigation")[0] : null;
    var place = placeTarget(Array.from(document.querySelectorAll("pre > code.sourceCode")));
    if (place && place.atTop) place = null;
    var moved = !!(nav && nav.type !== "navigate") && !place;
    // the browser restores its own offset around the load event: scrolls
    // until then are not the reader's, and the place is corrected after it
    var loaded = document.readyState === "complete", loadedAt = loaded ? Date.now() : 0, relands = 0;
    if (!loaded) window.addEventListener("load", function () {
      nextFrame(function () { nextFrame(function () { loaded = true; loadedAt = Date.now(); if (place) schedule(); }); });
    });
    var landed = null, since = Date.now();
    ["wheel", "touchmove", "keydown", "pointerdown"].forEach(function (ev) {
      window.addEventListener(ev, function () { moved = true; }, { passive: true, once: true });
    });
    window.addEventListener("hashchange", function () { moved = true; });
    function hashTarget() {
      if (!location.hash) return null;
      try { return document.getElementById(decodeURIComponent(location.hash.slice(1))); } catch (e) { return null; }
    }
    function landingEl() { return place ? place.el : hashTarget(); }
    window.addEventListener("scroll", function () {
      if (moved || landed === null || !loaded) return;
      var el = landingEl();
      if (el && Math.abs(el.getBoundingClientRect().top - landed) <= 48) return;
      // a returning page's saved place is re-landed when a scroll without
      // reader input moves it away just after load -- Chrome restores its
      // own offset just after the load event, over the correction -- but
      // only within 0.7 s of load and twice at most: later, a scroll without
      // input is find, focus or a screen reader, and is the reader's
      if (place && el && relands < 2 && loadedAt && Date.now() - loadedAt <= 700) { relands++; schedule(); }
      else moved = true;
    }, { passive: true });
    function refit() {
      charWidths = {};
      // (a refit keeps the reader's place; a page opened at a #section, or
      // returned to, re-lands below)
      if (moved || !landingEl()) preserving(fitSegments); else fitSegments();
      // (a page settles within seconds; after that the place is the reader's)
      if (moved || Date.now() - since > 10000) return;
      var el = landingEl();
      if (!el) return;
      if (place) window.scrollBy({ top: el.getBoundingClientRect().top - place.top, behavior: "instant" });
      else el.scrollIntoView({ block: "start", behavior: "instant" });
      landed = el.getBoundingClientRect().top;
    }
    function schedule() {
      if (queued) return;
      queued = true;
      nextFrame(function () { queued = false; refit(); });
    }
    refit();
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(schedule);
    // a font that arrives later, and a block whose width changes (shown from
    // a hidden tab or <details>, the column resized) refit too
    if (document.fonts && document.fonts.addEventListener) document.fonts.addEventListener("loadingdone", schedule);
    // On a long page a resize refits the blocks near the window at once and
    // the rest while the page is idle (a rotation refitted 3,000 lines in
    // one task); a newer resize drops batches still waiting.
    var near = new Set(), gen = 0, queuedR = false;
    if (typeof IntersectionObserver !== "undefined") {
      var nio = new IntersectionObserver(function (entries) {
        entries.forEach(function (e) {
          var c = e.target.querySelector(":scope > code.sourceCode");
          if (!c) return;
          if (e.isIntersecting) near.add(c); else near.delete(c);
        });
      }, { rootMargin: "100% 0px" });
      Array.from(document.querySelectorAll("pre > code.sourceCode")).forEach(function (c) { nio.observe(c.parentElement); });
    }
    function refitResized() {
      var landing = !moved && Date.now() - since <= 10000 && landingEl();
      if (!document.documentElement.hasAttribute("data-albers-long") || landing || typeof IntersectionObserver === "undefined") { refit(); return; }
      charWidths = {};
      var all = Array.from(document.querySelectorAll("pre > code.sourceCode[data-albers-seg]"));
      var now = all.filter(function (c) { return near.has(c); });
      var rest = all.filter(function (c) { return !near.has(c); });
      var g = ++gen;
      if (now.length) preserving(function () { fitSegments(now); });
      var batch = [], n = 0;
      function flush(b) {
        idle(function () { whenSettled(function () { if (g === gen) preserving(function () { fitSegments(b); }); }); });
      }
      rest.forEach(function (c) {
        batch.push(c); n += c.childElementCount;
        if (n >= CHUNK * 2) { flush(batch); batch = []; n = 0; }
      });
      if (batch.length) flush(batch);
    }
    function scheduleResized() {
      if (queuedR) return;
      queuedR = true;
      nextFrame(function () { queuedR = false; refitResized(); });
    }
    if (typeof ResizeObserver !== "undefined") {
      var widths = new WeakMap();
      var ro = new ResizeObserver(function (entries) {
        var changed = entries.some(function (e) {
          var w = Math.round(e.contentRect.width), was = widths.get(e.target);
          widths.set(e.target, w);
          return was !== undefined && was !== w;
        });
        if (changed) scheduleResized();
      });
      Array.from(document.querySelectorAll("pre > code.sourceCode")).forEach(function (c) { ro.observe(c.parentElement); });
    } else {
      window.addEventListener("resize", scheduleResized);
    }
    window.addEventListener("beforeprint", function () { printingCode = true; refit(); });
    window.addEventListener("afterprint", function () { printingCode = false; refit(); });
  }

  // Each source line hangs from its own indentation when it wraps.
  // Restructure a code element while it is out of the document: moving
  // thousands of nodes inside a live block re-checks the page's :has()
  // selectors on every move (seconds on a long listing).
  function detached(code, fn) {
    var parent = code.parentNode, next = code.nextSibling;
    parent.removeChild(code);
    try { fn(); } finally { parent.insertBefore(code, next); }
  }

  // A very long block is prepared a screenful at a time: its first lines
  // now, the rest in later tasks (one 1,500-line listing held the page for
  // over half a second). finishCode completes one at once, before a jump
  // past it or printing.
  var CHUNK = 300;
  function segmentLines(code, lines) {
    detached(code, function () {
      lines.forEach(function (line) {
        var m = /^( +)/.exec(line.textContent);
        if (m) { line.style.setProperty("--ind", m[1].length + "ch"); line.setAttribute("data-ind", m[1].length); }
        segmentLine(line);
      });
    });
  }
  function segmentCode(code) {
    if (code.hasAttribute("data-albers-seg")) return;
    code.setAttribute("data-albers-seg", "");
    var lines = Array.from(code.querySelectorAll(":scope > span:not(.albers-out-block)"));
    if (lines.length <= CHUNK * 1.5) { segmentLines(code, lines); return; }
    segmentLines(code, lines.slice(0, CHUNK));
    code._albersRest = lines.slice(CHUNK);
    // the next chunk when its first line comes within reach (preparing
    // every chunk in idle time blocked the page for over a second)
    if (typeof IntersectionObserver === "undefined") {
      idle(function step() {
        var rest = code._albersRest;
        if (!rest || !rest.length) return;
        var batch = rest.splice(0, CHUNK);
        preserving(function () { segmentLines(code, batch); fitSegments(code, batch); });
        if (rest.length) idle(step);
      });
      return;
    }
    // every chunk's first line is watched, so a jump into the middle of the
    // block (find, a line link) prepares the chunk it lands in -- and the
    // ones before it, so the landing does not move
    // (every 50th line is a sentinel for its chunk: a chunk is taller than
    // the window, so a jump into its middle reached neither end)
    var chunks = [], owner = new Map();
    for (var ci = 0; ci < code._albersRest.length; ci += CHUNK) chunks.push(code._albersRest.slice(ci, ci + CHUNK));
    chunks.forEach(function (c, k) { for (var j = 0; j < c.length; j += 50) owner.set(c[j], k); });
    var sio = new IntersectionObserver(function (entries) {
      var hit = -1, hitEl = null;
      entries.forEach(function (e) {
        if (e.isIntersecting && owner.has(e.target) && owner.get(e.target) >= hit) { hit = owner.get(e.target); hitEl = e.target; }
      });
      if (hit < 0) return;
      var run = function () {
        var todo = [];
        for (var k = 0; k <= hit; k++) if (chunks[k]) {
          todo = todo.concat(chunks[k]);
          for (var j = 0; j < chunks[k].length; j += 50) sio.unobserve(chunks[k][j]);
          chunks[k] = null;
        }
        var left = new Set(code._albersRest || []);
        todo = todo.filter(function (ln) { return left.has(ln); });
        if (!todo.length) return;
        var picked = new Set(todo);
        code._albersRest = code._albersRest.filter(function (ln) { return !picked.has(ln); });
        preserving(function () { segmentLines(code, todo); fitSegments(code, todo); });
      };
      whenSettled(run, hitEl, function () { sio.unobserve(hitEl); sio.observe(hitEl); });
    }, { rootMargin: "1200px 0px" });
    owner.forEach(function (k, ln) { sio.observe(ln); });
  }

  // Work that can wait runs when the page is idle (not before the first
  // paint, and not in the middle of a scroll), one piece per slice.
  var idleQ = [], idleOn = false;
  // (a frame where frames run; a short timer elsewhere)
  function nextFrame(fn) {
    return window.requestAnimationFrame ? window.requestAnimationFrame(fn) : setTimeout(fn, 16);
  }
  function idle(fn) {
    idleQ.push(fn);
    if (!idleOn) { idleOn = true; pump(); }
  }
  // (without requestIdleCallback -- Safari -- a timer slice of about 8 ms;
  // a job that throws does not stop the queue; nothing runs while a mouse
  // button is down, so a drag-selection is not restructured under it)
  var pointerDown = false, pointerTimer = null;
  // (a release that never arrives -- a context menu -- counts after 3 s)
  document.addEventListener("pointerdown", function () {
    pointerDown = true;
    clearTimeout(pointerTimer);
    pointerTimer = setTimeout(function () { pointerDown = false; if (idleQ.length && !idleOn) { idleOn = true; pump(); } }, 3000);
  }, true);
  document.addEventListener("pointerup", function () { pointerDown = false; }, true);
  document.addEventListener("pointercancel", function () { pointerDown = false; }, true);
  // A scroll the reader did not drive (a smooth scroll to a link target,
  // find, scrollIntoView) is not interrupted: preparing code with its place
  // correction cancelled it, stopping far short. Work that would move the
  // page waits until that scroll ends.
  // (the reader's own scrolling: the wheel, touch, the scrolling keys, and a
  // press that is not on a link -- a link's press starts the very smooth
  // scroll that must not be interrupted, and counted as input it let
  // preparation cancel a contents click mid-flight)
  var lastInput = 0, autoScrolling = false, autoEnd = null, afterScroll = [];
  var SCROLL_KEYS = /^(ArrowUp|ArrowDown|PageUp|PageDown|Home|End| |Spacebar)$/;
  function readerInput(e) {
    if (e.type === "keydown" && !SCROLL_KEYS.test(e.key)) return;
    if (e.type === "pointerdown" && e.target && e.target.closest && e.target.closest("a[href]")) return;
    lastInput = Date.now();
  }
  ["wheel", "touchmove", "keydown", "pointerdown"].forEach(function (ev) {
    window.addEventListener(ev, readerInput, { passive: true, capture: true });
  });
  function endAutoScroll() {
    autoScrolling = false;
    var todo = afterScroll; afterScroll = [];
    // only what is still near the window runs; the rest is handed back to
    // its observer (a smooth scroll passes every block on its way)
    var reach = window.innerHeight + 1200;
    var near = todo.map(function (t) {
      if (!t.el) return true;
      var r = t.el.getBoundingClientRect();
      return r.bottom > -1200 && r.top < reach;
    });
    todo.forEach(function (t, i) {
      try { if (near[i]) t.fn(); else if (t.requeue) t.requeue(); } catch (e) {}
    });
  }
  window.addEventListener("scroll", function () {
    if (Date.now() - lastInput < 250) return;
    autoScrolling = true;
    clearTimeout(autoEnd);
    autoEnd = setTimeout(endAutoScroll, 160);
  }, { passive: true });
  if ("onscrollend" in window) window.addEventListener("scrollend", function () {
    if (autoScrolling) { clearTimeout(autoEnd); endAutoScroll(); }
  });
  // run fn now, or when a drag or an automatic scroll is over (el: what fn
  // prepares, checked for nearness then; requeue: hands it back otherwise)
  function whenSettled(fn, el, requeue) {
    if (autoScrolling) afterScroll.push({ fn: fn, el: el, requeue: requeue });
    else if (pointerDown) idle(fn);
    else fn();
  }

  function pump() {
    var ric = window.requestIdleCallback || function (cb) {
      return setTimeout(function () {
        var end = Date.now() + 8;
        cb({ timeRemaining: function () { return Math.max(0, end - Date.now()); } });
      }, 60);
    };
    ric(function (deadline) {
      if (!pointerDown) {
        do {
          try { idleQ.shift()(); } catch (e) {}
        } while (idleQ.length && deadline.timeRemaining() > 8);
      }
      if (idleQ.length) pump(); else idleOn = false;
    }, { timeout: 1500 });
  }

  // Changing code the reader may be looking at keeps their place (the element
  // at the top of the window stays where it was) and their selection (as
  // character offsets in the block, since text nodes are split). A find-bar
  // match cannot be kept, so blocks are prepared soon after load.
  function preserving(fn) {
    var anchor = null, top = 0, sel = null;
    if (window.scrollY > 0 && document.elementFromPoint) {
      anchor = document.elementFromPoint(Math.round(window.innerWidth / 2), Math.round(Math.min(window.innerHeight / 3, 160)));
      if (anchor) top = anchor.getBoundingClientRect().top;
    }
    // each end of the selection on its own (it may span blocks), and its
    // direction (anchor and focus, not start and end)
    var s = window.getSelection && window.getSelection();
    if (s && s.rangeCount && !s.isCollapsed && s.anchorNode && s.focusNode) {
      sel = { a: endpoint(s.anchorNode, s.anchorOffset), f: endpoint(s.focusNode, s.focusOffset) };
    }
    fn();
    if (sel && (sel.a.code || sel.f.code)) {
      var a = restorePoint(sel.a), b = restorePoint(sel.f);
      if (a && b) { try { s.setBaseAndExtent(a[0], a[1], b[0], b[1]); } catch (e) {} }
    }
    if (anchor && anchor.isConnected) {
      var d = anchor.getBoundingClientRect().top - top;
      if (Math.abs(d) > 0.5) window.scrollBy({ top: d, behavior: "instant" });
    }
  }
  function endpoint(node, off) {
    var el = node.nodeType === 1 ? node : node.parentNode;
    var code = el && el.closest && el.closest("pre > code");
    return code ? { code: code, n: offsetIn(code, node, off) } : { node: node, off: off };
  }
  function restorePoint(p) {
    return p.code ? pointAt(p.code, p.n) : (p.node.isConnected ? [p.node, p.off] : null);
  }
  function offsetIn(root, node, off) {
    var r = document.createRange();
    r.setStart(root, 0);
    r.setEnd(node, off);
    return r.toString().length;
  }
  function pointAt(root, n) {
    var tw = document.createTreeWalker(root, NodeFilter.SHOW_TEXT), t;
    while ((t = tw.nextNode())) {
      if (n <= t.nodeValue.length) return [t, n];
      n -= t.nodeValue.length;
    }
    return null;
  }
  function finishCode(code) {
    segmentCode(code);
    var rest = code._albersRest;
    if (rest && rest.length) segmentLines(code, rest.splice(0, rest.length));
  }


  // Blocks near the top are prepared now (so the first paint is final); the
  // rest as they come within reach of the viewport, and all before printing.
  // A page opened at a #section prepares every block before that section
  // first, so the section is where the reader lands (lazy blocks above it
  // once moved it). Positions are read in one pass, before any block changes.
  // The reader's place: the first code line, or prose block (heading,
  // paragraph, list item, figure...), found down the reading line. A code
  // line is kept as (block, line); anything else as its index among the
  // page's blocks, which preparing code never changes. (A single point
  // fell in the gap between blocks and hit <main>: the place was lost.)
  var PLACE_SEL = ":is(h1, h2, h3, h4, h5, h6, p, li, dt, dd, figure, table, blockquote, pre, .callout)";
  function placeBlocks() {
    // (pkgdown's <main>; a vignette has none: its body)
    var scope = document.querySelector("main") || document.body;
    return Array.prototype.slice.call(scope.querySelectorAll(PLACE_SEL));
  }
  // a code line's text with the lines around it (repeated lines differ in
  // their neighbours)
  function lineContext(line) {
    var prev = line.previousElementSibling, next = line.nextElementSibling;
    return [prev, line, next].map(function (l) { return l ? placeText(l).slice(0, 40) : ""; }).join("|");
  }
  function placeText(el) { return (el.textContent || "").replace(/\s+/g, " ").trim().slice(0, 60); }
  function savePlace(codes) {
    try {
      var rec = null, x = [0.5, 0.3, 0.7], y0 = Math.round(Math.min(window.innerHeight / 3, 160));
      for (var dy = 0; !rec && dy <= 240; dy += 20) {
        for (var i = 0; !rec && i < x.length; i++) {
          var el = document.elementFromPoint(Math.round(window.innerWidth * x[i]), y0 + dy);
          if (!el || !el.closest) continue;
          // (not the navigation, a sticky contents list or a header)
          if (el.closest("nav, #TOC, #toc, aside, header, .navbar")) continue;
          var line = el.closest("pre > code.sourceCode > span");
          if (line) {
            var code = line.parentElement;
            rec = { c: codes.indexOf(code), l: Array.prototype.indexOf.call(code.children, line), top: line.getBoundingClientRect().top, t: lineContext(line) };
            continue;
          }
          var block = el.closest(PLACE_SEL);
          if (block) {
            var k = placeBlocks().indexOf(block);
            if (k >= 0) rec = { b: k, top: block.getBoundingClientRect().top, t: placeText(block) };
          }
        }
      }
      // (at the top: nothing to correct, and nothing to prepare on return)
      if (window.scrollY <= 0) rec = { atTop: true };
      // kept in the history entry itself (storage may be blocked) and in
      // session storage
      // (only into a state that is null or a plain object -- another script's
      // string or array is left alone -- only when the place changed, and
      // not into a very large one)
      try {
        var hs = history.state;
        var plain = hs === null || (typeof hs === "object" && !Array.isArray(hs) && Object.getPrototypeOf(hs) === Object.prototype);
        if (plain) {
          var before = hs && hs.albersPlace ? JSON.stringify(hs.albersPlace) : "null";
          if (before !== JSON.stringify(rec || null) && (!hs || JSON.stringify(hs).length < 65536)) {
            history.replaceState(Object.assign({}, hs || {}, { albersPlace: rec || null }), "");
          }
        }
      } catch (e3) {}
      try {
        if (rec) sessionStorage.setItem(PLACE_KEY, JSON.stringify(rec));
        else sessionStorage.removeItem(PLACE_KEY);
      } catch (e4) {}
    } catch (e) {}
  }
  function placeTarget(codes) {
    if (!savedPlace) return null;
    var el = null;
    if (savedPlace.atTop) return { atTop: true };
    if (savedPlace.b !== undefined) {
      // (by index, checked against its text: a page rebuilt between visits
      // shifts indices, so the nearest block with the same text is taken)
      var blocks = placeBlocks();
      el = blocks[savedPlace.b] || null;
      if (savedPlace.t && (!el || placeText(el) !== savedPlace.t)) {
        el = null;
        for (var d = 1; d < blocks.length && !el; d++) {
          [savedPlace.b - d, savedPlace.b + d].forEach(function (k) {
            if (!el && blocks[k] && placeText(blocks[k]) === savedPlace.t) el = blocks[k];
          });
        }
      }
    }
    else if (savedPlace.id) el = document.getElementById(savedPlace.id);
    else if (savedPlace.c !== undefined) {
      // (a code line by block and line, checked against its text; if the
      // page was rebuilt, the nearest line with the same text)
      el = codes[savedPlace.c] ? codes[savedPlace.c].children[savedPlace.l] || null : null;
      if (savedPlace.t && (!el || lineContext(el) !== savedPlace.t)) {
        // (the line with its neighbours, so a repeated line is not taken for
        // another; if nothing matches -- the line was edited -- its old place
        // by index is the best guess left)
        var byIndex = el;
        el = null;
        var lines = [];
        codes.forEach(function (c) { Array.prototype.push.apply(lines, c.children); });
        var near = codes[savedPlace.c] && codes[savedPlace.c].children[0] ? lines.indexOf(codes[savedPlace.c].children[0]) + savedPlace.l : 0;
        for (var dd = 0; dd < lines.length && !el; dd++) {
          [near - dd, near + dd].forEach(function (k) {
            if (!el && lines[k] && lineContext(lines[k]) === savedPlace.t) el = lines[k];
          });
        }
        el = el || byIndex;
      }
    }
    return el ? { el: el, top: savedPlace.top || 0 } : null;
  }

  function initLineIndents() {
    var codes = Array.from(document.querySelectorAll("pre > code.sourceCode"));
    var total = codes.reduce(function (n, c) { return n + c.childElementCount; }, 0);
    // a typical vignette is prepared at once
    // (on a reload or Back/Forward the browser restores a pixel offset saved
    // with every block prepared: prepare them all, or the place is lost)
    var nav = performance.getEntriesByType ? performance.getEntriesByType("navigation")[0] : null;
    var restoring = !!(nav && nav.type !== "navigate");
    if (total > 600) {
      document.documentElement.setAttribute("data-albers-long", "");
      window.addEventListener("pagehide", function () { savePlace(codes); });
      // (and a moment after scrolling stops: a history entry cannot always
      // be written while the page unloads)
      // (not while a return is still being corrected: the first 2 s)
      var placeTimer = null, placeFrom = Date.now() + 2000;
      window.addEventListener("scroll", function () {
        clearTimeout(placeTimer);
        placeTimer = setTimeout(function () { if (Date.now() >= placeFrom) savePlace(codes); }, 400);
      }, { passive: true });
    }
    // (a restore without a saved place: prepare everything, so the
    // browser's saved offset still means what it meant)
    var place = restoring && total > 600 ? placeTarget(codes) : null;
    // (a saved place that is no longer on the page -- edited between visits
    // -- is not a reason to prepare everything: that cost two seconds)
    if (typeof IntersectionObserver === "undefined" || total <= 600 || (restoring && !place && !savedPlace)) { codes.forEach(finishCode); return; }
    if (place && place.atTop) place = null;
    var target = place ? place.el : null;
    try { if (!target && location.hash) target = document.getElementById(decodeURIComponent(location.hash.slice(1))); } catch (e) {}
    if (place) {
      // the block holding the saved line is prepared whole too
      var holder = place.el.closest && place.el.closest("pre > code.sourceCode");
      if (holder) finishCode(holder);
    }
    // the first blocks, by line count (a screenful or two of code), are
    // prepared now; reading their positions would cost a layout of the page
    var budget = 400;
    var later = [];
    preserving(function () {
      codes.forEach(function (code) {
        var before = target && (target.compareDocumentPosition(code) & Node.DOCUMENT_POSITION_PRECEDING);
        if (before) finishCode(code);
        else if (budget > 0) { budget -= code.childElementCount; segmentCode(code); }
        else later.push(code);
      });
    });
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (!e.isIntersecting) return;
        io.unobserve(e.target);
        var code = e.target.querySelector(":scope > code.sourceCode");
        if (!code || code.hasAttribute("data-albers-seg")) return;
        // (not under a drag-selection or an automatic scroll in progress)
        var pre = e.target;
        whenSettled(function () {
          if (!code.hasAttribute("data-albers-seg")) preserving(function () { segmentCode(code); fitSegments(code); });
        }, pre, function () { io.observe(pre); });
      });
    }, { rootMargin: "1200px 0px" });
    // (only as they come near: preparing every block in idle time blocked
    // the page for seconds on a long article)
    later.forEach(function (c) { io.observe(c.parentElement); });
    // (a jump within the page no longer prepares the blocks before its
    // target first -- that cost a second per contents click on a long page;
    // blocks prepared around the landing keep it in place instead)
    window.addEventListener("beforeprint", function () {
      codes.forEach(finishCode); fitSegments(codes);
    });
  }

  function initWideBlocks() {
    Array.from(document.querySelectorAll("pre.wide")).forEach(function (pre) {
      if (pre.parentElement && pre.parentElement.matches("div.sourceCode")) pre.parentElement.classList.add("wide");
    });
  }

  // Both images stay in the figure; CSS shows the dark one only on screen in
  // dark mode, so printing always uses the light figure.
  function initFigureTwins() {
    Array.from(document.querySelectorAll("img.albers-dark-twin")).forEach(function (twin) {
      var holder = twin.parentElement && twin.parentElement.tagName === "P" && twin.parentElement.children.length === 1
        ? twin.parentElement : twin;
      var prev = holder.previousElementSibling;
      var img = prev && (prev.tagName === "IMG" ? prev : prev.querySelector("img:not(.albers-dark-twin)"));
      if (img && !img.classList.contains("albers-has-twin")) {
        img.classList.add("albers-has-twin");
        twin.removeAttribute("hidden");
        twin.setAttribute("loading", "lazy");
        twin.setAttribute("alt", img.getAttribute("alt") || "");
        twin.removeAttribute("aria-hidden");
        // the same size and placement as the light image (out.width, fig.align)
        ["width", "height"].forEach(function (a) {
          if (img.hasAttribute(a)) twin.setAttribute(a, img.getAttribute(a));
        });
        if (img.style.margin) twin.style.margin = img.style.margin;
        if (img.style.width) twin.style.width = img.style.width;
        // shown the way its light image is (side-by-side figures are inline)
        twin.style.setProperty("--albers-twin-display", getComputedStyle(img).display === "inline" ? "inline-block" : getComputedStyle(img).display);
        img.parentNode.insertBefore(twin, img.nextSibling);
      }
      if (holder !== twin && holder.parentNode) holder.parentNode.removeChild(holder);
    });
  }

  // Phone twins: albers_vignette() also draws each plot at phone width and
  // puts it, hidden, after the figure (the dark twin's in the paragraph after
  // it). While a figure is shown narrower than data-albers-below (CSS px), it
  // shows the phone drawing, whose text keeps its size; print and the zoom
  // use the full drawing.
  function initPhoneTwins() {
    var swaps = [];
    Array.from(document.querySelectorAll("img.albers-phone-twin")).forEach(function (ph) {
      var holder = ph.parentElement && ph.parentElement.tagName === "P" && ph.parentElement.children.length === 1
        ? ph.parentElement : ph;
      var img;
      if (ph.classList.contains("albers-phone-dark")) {
        img = holder.previousElementSibling;
        if (img && !img.matches("img.albers-dark-twin")) img = img.querySelector("img.albers-dark-twin");
      } else {
        img = ph.previousElementSibling;
        while (img && img.matches("img.albers-dark-twin, img.albers-phone-twin")) img = img.previousElementSibling;
      }
      var below = parseFloat(ph.getAttribute("data-albers-below"));
      if (img && img.tagName === "IMG" && below > 0) {
        img.albersFull = img.getAttribute("src");
        swaps.push({ img: img, phone: ph.getAttribute("src"), below: below });
      }
      if (holder.parentNode) holder.parentNode.removeChild(holder);
    });
    if (!swaps.length) return;
    var printing = false, queued = false;
    function fit() {
      queued = false;
      swaps.forEach(function (s) {
        // the light image and its dark twin share one width; one is hidden
        var w = 0;
        Array.from(s.img.parentElement.children).forEach(function (el) {
          if (el.tagName === "IMG") w = Math.max(w, el.getBoundingClientRect().width);
        });
        // (only on a phone-width screen: a figure shown narrow on a desktop
        // -- out.width 50%, side by side -- keeps its own drawing and shape)
        var phoneScreen = !!(window.matchMedia && window.matchMedia("(max-width: 620px)").matches);
        var src = !printing && phoneScreen && w > 0 && w < s.below ? s.phone : s.img.albersFull;
        if (s.img.getAttribute("src") !== src) s.img.setAttribute("src", src);
      });
    }
    fit();
    function queue() { if (!queued) { queued = true; nextFrame(fit); } }
    window.addEventListener("resize", queue);
    // (a figure first shown later -- in <details>, a tab -- changes size)
    if (typeof ResizeObserver !== "undefined") {
      var fro = new ResizeObserver(queue);
      swaps.forEach(function (s) { fro.observe(s.img); });
    }
    window.addEventListener("beforeprint", function () { printing = true; fit(); });
    window.addEventListener("afterprint", function () { printing = false; fit(); });
  }

  /* pkgdown article headers carry the same title plate as vignettes */
  // pkgdown pages: the same names and states the vignette script provides.
  function initSiteA11y() {
    // section permalinks are named after their section, not "anchor"
    Array.from(document.querySelectorAll("a.anchor[href^='#']")).forEach(function (a) {
      var h = a.closest("h1, h2, h3, h4, h5, h6");
      if (h) a.setAttribute("aria-label", "Link to section: " + (h.textContent || "").trim());
    });
    // footnote references: "Footnote 1", and Escape closes the note
    Array.from(document.querySelectorAll("a.footnote-ref")).forEach(function (a) {
      if (!a.hasAttribute("aria-label")) a.setAttribute("aria-label", "Footnote " + (a.textContent || "").trim());
    });
    document.addEventListener("keydown", function (e) {
      if (e.key !== "Escape" || !window.bootstrap || !bootstrap.Popover) return;
      Array.from(document.querySelectorAll("[data-bs-toggle='popover']")).forEach(function (el) {
        var pop = bootstrap.Popover.getInstance(el);
        if (pop) pop.hide();
      });
    });
    // contents: Bootstrap's scrollspy marks the entry and its parent active;
    // the deepest one is the current location
    var toc = document.getElementById("toc");
    if (!toc) return;
    function sync() {
      var active = Array.from(toc.querySelectorAll("a.nav-link.active"));
      var deepest = active.filter(function (a) { return !a.parentElement.querySelector(":scope ul a.nav-link.active"); }).pop();
      Array.from(toc.querySelectorAll("a[aria-current]")).forEach(function (a) { if (a !== deepest) a.removeAttribute("aria-current"); });
      if (deepest && deepest.getAttribute("aria-current") !== "location") deepest.setAttribute("aria-current", "location");
    }
    if (typeof MutationObserver === "undefined") return;
    new MutationObserver(sync).observe(toc, { subtree: true, attributes: true, attributeFilter: ["class"] });
    sync();
    // (the location mark is albers.js's: without it, CSS mutes nothing)
    toc.classList.add("albers-toc-sync");
  }

  function initSitePlate() {
    var header = document.querySelector(".template-article main#main > .page-header, .template-home main#main .page-header");
    if (!header || header.querySelector(".albers-plate-mark") || document.body.classList.contains("albers-no-plate")) return;
    header.classList.add("albers-titleblock", "albers-titleblock--site");
    var plate = document.createElement("div");
    plate.className = "albers-plate-mark";
    plate.setAttribute("aria-hidden", "true");
    header.appendChild(plate);
  }

  /* ---------------------------------------------------------------------
   * Contents: landmark, scrollspy, collapsible on narrow screens
   * ------------------------------------------------------------------- */
  function initToc() {
    var toc = document.getElementById("TOC");
    if (!toc) return;
    toc.setAttribute("role", "navigation");
    toc.setAttribute("aria-label", "Contents");

    var list = toc.querySelector(":scope > ul");
    if (list && !toc.querySelector(".albers-toc-toggle")) {
      var btn = document.createElement("button");
      btn.type = "button";
      btn.className = "albers-toc-toggle";
      btn.textContent = "Contents";
      btn.setAttribute("aria-expanded", "false");
      if (!list.id) list.id = "albers-toc-list";
      btn.setAttribute("aria-controls", list.id);
      btn.addEventListener("click", function () {
        var open = toc.classList.toggle("is-open");
        btn.setAttribute("aria-expanded", String(open));
      });
      toc.insertBefore(btn, list);
      toc.classList.add("has-toggle");
      function closePanel() {
        if (!toc.classList.contains("is-open")) return;
        toc.classList.remove("is-open");
        btn.setAttribute("aria-expanded", "false");
      }
      toc.addEventListener("click", function (e) {
        if (e.target.closest && e.target.closest("a") && window.matchMedia && window.matchMedia("(max-width: 1099.98px)").matches) closePanel();
      });
      // the open panel closes when focus or a tap leaves it
      toc.addEventListener("focusout", function (e) {
        if (e.relatedTarget && !toc.contains(e.relatedTarget)) closePanel();
      });
      document.addEventListener("pointerdown", function (e) {
        if (!toc.contains(e.target)) closePanel();
      });
    }

    var links = Array.from(toc.querySelectorAll("a[href^='#']"));
    if (!links.length || typeof IntersectionObserver === "undefined") return;
    var byId = {};
    var targets = [];
    links.forEach(function (a) {
      var id = decodeURIComponent(a.getAttribute("href").slice(1));
      var el = document.getElementById(id);
      if (!el) return;
      byId[id] = a;
      targets.push(el);
    });
    if (!targets.length) return;

    var visible = {};
    function setActive(id) {
      links.forEach(function (a) {
        a.classList.remove("is-active");
        a.removeAttribute("aria-current");
      });
      var a = byId[id];
      if (a) {
        a.classList.add("is-active");
        a.setAttribute("aria-current", "location");
      }
    }
    // A clicked entry, or the section named in the URL (on load, Back and
    // Forward), stays current while its heading is on screen, even when the
    // scroll position alone would name another (short sections at the end).
    var clicked = null, clickedAt = 0;
    function pin(id) {
      if (!byId[id]) return;
      clicked = id;
      clickedAt = Date.now();
      setActive(id);
    }
    function hashId() {
      try { return decodeURIComponent(location.hash.slice(1)); } catch (e) { return ""; }
    }
    links.forEach(function (a) {
      a.addEventListener("click", function () { pin(decodeURIComponent(a.getAttribute("href").slice(1))); });
    });
    window.addEventListener("hashchange", function () { pin(hashId()); });
    window.addEventListener("popstate", function () { if (!location.hash) { clicked = null; update(); } });
    if (location.hash) pin(hashId());
    // the reader scrolling on their own releases the pin
    function release() { clicked = null; }
    window.addEventListener("wheel", release, { passive: true });
    window.addEventListener("touchmove", release, { passive: true });
    window.addEventListener("keydown", function (e) {
      if (/^(ArrowUp|ArrowDown|PageUp|PageDown|Home|End| )$/.test(e.key)) release();
    });
    // Headings, not whole sections, decide which entry is current: the last
    // heading that has scrolled past the top third of the viewport.
    var heads = targets.map(function (el) {
      return el.matches("h1,h2,h3,h4") ? el : (el.querySelector(":scope > h1, :scope > h2, :scope > h3, :scope > h4") || el);
    });
    function update() {
      if (clicked) {
        var ph = heads[targets.map(function (t) { return t.id; }).indexOf(clicked)];
        var top = ph ? ph.getBoundingClientRect().top : -1;
        if (Date.now() - clickedAt < 1000 || (top >= -4 && top < window.innerHeight * 0.85)) { setActive(clicked); return; }
        clicked = null;
      }
      var line = window.innerHeight * 0.3;
      var current = null;
      for (var i = 0; i < heads.length; i++) {
        if (heads[i].getBoundingClientRect().top <= line) current = targets[i].id;
      }
      if (window.innerHeight + window.scrollY >= document.documentElement.scrollHeight - 4) {
        current = targets[targets.length - 1].id;
      }
      setActive(current);
    }
    var ticking = false;
    window.addEventListener("scroll", function () {
      if (ticking) return;
      ticking = true;
      nextFrame(function () { ticking = false; update(); });
    }, { passive: true });
    update();
  }

  /* ---------------------------------------------------------------------
   * Heading anchors. html_vignette puts ids on the section <div>, not on the
   * heading, so fall back to the enclosing section's id.
   * ------------------------------------------------------------------- */
  function initAnchors() {
    Array.from(document.querySelectorAll("h2, h3")).forEach(function (h) {
      if (h.querySelector("a.anchor, a.albers-heading-link") || h.closest("aside, nav, #TOC, .albers-titleblock")) return;
      var id = h.id;
      if (!id) {
        var sec = h.parentElement;
        if (sec && sec.id && sec.firstElementChild === h) id = sec.id;
      }
      if (!id) return;
      // The heading text itself becomes the permalink: operable by keyboard,
      // and the heading's accessible name stays "Overview" (the "#" cue is
      // CSS-generated, so it is not read out).
      if (h.querySelector("a")) return; // headings that already contain links
      var a = document.createElement("a");
      a.href = "#" + id;
      a.className = "albers-heading-link";
      while (h.firstChild) a.appendChild(h.firstChild);
      if (h.hasAttribute("aria-label")) a.setAttribute("aria-label", h.getAttribute("aria-label"));
      h.appendChild(a);
    });
  }

  /* ---------------------------------------------------------------------
   * Output lines. With collapse = TRUE knitr emits output as "#>" comments;
   * tag those lines so they read as results (and warnings / errors as such),
   * not as faint italic comments.
   * ------------------------------------------------------------------- */
  // albers_vignette() marks condition lines after the prompt, whatever the
  // chunk's comment prefix: U+2063 message, U+2064 warning, U+2062 error
  var MARKS = { "\u2063": "msg", "\u2064": "warn", "\u2062": "err" };
  var MARK_RE = /[\u2062\u2063\u2064]/;
  function markedKind(text) {
    var m = MARK_RE.exec(text);
    return m && m.index <= 24 ? MARKS[m[0]] : null;
  }
  function lineKind(text) {
    var t = text.replace(/^\s+/, "");
    var marked = markedKind(t);
    if (marked) return marked;
    if (t.indexOf("#>") !== 0) return null;
    var body = t.slice(2).replace(/^\s+/, "");
    if (/^(Warning|Warning message|Warning in)\b/.test(body)) return "warn";
    if (/^(Error|Error in)\b/.test(body) || /^!\s/.test(body)) return "err";
    return "out";
  }

  function wrapPrompt(line) {
    var walker = document.createTreeWalker(line, NodeFilter.SHOW_TEXT);
    var node;
    while ((node = walker.nextNode())) {
      var i = node.nodeValue.indexOf("#>");
      if (i === -1) continue;
      var after = node.splitText(i);
      var rest = after.splitText(2);
      var span = document.createElement("span");
      span.className = "albers-prompt";
      span.textContent = "#>";
      after.parentNode.replaceChild(span, after);
      if (rest.nodeValue.charAt(0) === " ") {
        // keep the separating space with the (unselectable) prompt
        rest.nodeValue = rest.nodeValue.slice(1);
        span.textContent = "#> ";
      }
      return;
    }
  }

  function stripMark(line) {
    var walker = document.createTreeWalker(line, NodeFilter.SHOW_TEXT);
    var node;
    while ((node = walker.nextNode())) {
      if (MARK_RE.test(node.nodeValue)) node.nodeValue = node.nodeValue.replace(/[\u2062\u2063\u2064]/g, "");
    }
  }

  // Output that knitr did not collapse into the code block (collapse = FALSE,
  // echo = FALSE) is a plain <pre>: a message there wraps like prose. Then no
  // U+2063 mark is left anywhere to be copied, found or printed.
  function initPlainMessages() {
    // an echo-free chunk is a plain <pre> of "#>" lines: it is output, and
    // reads as output (the ink and pitch of an output band)
    Array.from(document.querySelectorAll("pre:not(.sourceCode) > code:not(.sourceCode)")).forEach(function (code) {
      var lines = code.textContent.split("\n").filter(function (l) { return /\S/.test(l); });
      // ("#>" by default; "##" is knitr's own default comment)
      if (lines.length && lines.every(function (l) { return /^(#>|##)/.test(l); })) code.parentElement.classList.add("albers-plain-out");
    });
    Array.from(document.querySelectorAll("pre > code")).forEach(function (code) {
      if (!MARK_RE.test(code.textContent)) return;
      if (!code.querySelector(":scope > span")) {
        var lines = code.textContent.split("\n").filter(function (l) { return l.length; });
        var kinds = lines.map(markedKind);
        var uniform = kinds.every(function (k) { return k && k === kinds[0]; });
        if (uniform) {
          code.parentElement.classList.add("albers-msg-pre", "is-" + kinds[0]);
        } else {
          // a condition after printed output in one block: only its lines wrap
          var parts = code.textContent.split("\n");
          code.textContent = "";
          parts.forEach(function (l, i) {
            if (i) code.appendChild(document.createTextNode("\n"));
            var k = markedKind(l);
            if (k) {
              var span = document.createElement("span");
              span.className = "albers-msg-line is-" + k;
              span.textContent = l;
              code.appendChild(span);
            } else if (l.length) {
              code.appendChild(document.createTextNode(l));
            }
          });
        }
      }
      var walker = document.createTreeWalker(code, NodeFilter.SHOW_TEXT);
      var n;
      while ((n = walker.nextNode())) {
        if (MARK_RE.test(n.nodeValue)) n.nodeValue = n.nodeValue.replace(/[\u2062\u2063\u2064]/g, "");
      }
    });
    // pandoc's per-line anchors are empty links: out of the accessibility tree
    Array.from(document.querySelectorAll("pre code a:empty[href^='#']")).forEach(function (a) {
      a.setAttribute("aria-hidden", "true");
      a.setAttribute("tabindex", "-1");
    });
  }

  // the format's post-processor turns condition marks into empty
  // <span class="albers-mk" data-k="..."> elements; put the characters back
  // for the passes below (which remove them again)
  function unpackMarks() {
    var ch = { msg: "\u2063", warn: "\u2064", err: "\u2062" };
    Array.from(document.querySelectorAll("span.albers-mk")).forEach(function (mk) {
      var c = ch[mk.getAttribute("data-k")];
      if (c) mk.parentNode.replaceChild(document.createTextNode(c), mk);
      else mk.parentNode.removeChild(mk);
    });
    // merge the split text nodes, so a line's text reads as one
    Array.from(document.querySelectorAll("pre code")).forEach(function (c) { c.normalize(); });
  }

  function initOutputLines() {
    unpackMarks();
    // on a page rendered by albers_vignette() every condition line is marked,
    // so an unmarked "#> Error ..." is printed output, not an error
    var pageMarked = document.documentElement.hasAttribute("data-albers-marks");
    Array.from(document.querySelectorAll("pre code")).forEach(function (code) {
      if (code.querySelector(".albers-out-block")) return;
      var lines = Array.from(code.children).filter(function (c) { return c.tagName === "SPAN"; });
      if (!lines.length) return;

      // 1. classify lines. Marked lines (albers_vignette) say exactly what
      // they are; without marks, text after a warning/error line is taken as
      // its wrapped continuation unless it starts a printed vector.
      var isMarked = pageMarked || MARK_RE.test(code.textContent);
      var kinds = [];
      var prevKind = null;
      lines.forEach(function (line) {
        var kind = lineKind(line.textContent);
        if (isMarked) {
          kind = markedKind(line.textContent) || (kind ? "out" : null);
        } else if (kind === "out" && (prevKind === "warn" || prevKind === "err") && !/^\s*#>\s*\[\d+\]/.test(line.textContent)) {
          kind = prevKind;
        }
        kinds.push(kind);
        prevKind = kind;
      });

      // 2. gather runs of output lines into one block each
      var i = 0;
      while (i < lines.length) {
        if (!kinds[i]) { i++; continue; }
        var j = i;
        while (j + 1 < lines.length && kinds[j + 1]) j++;
        var block = document.createElement("span");
        block.className = "albers-out-block";
        var runKinds = kinds.slice(i, j + 1);
        // the block takes a condition's band only when every line is that
        // condition; in a mixed run each line carries its own kind
        var uniform = runKinds.every(function (k) { return k === runKinds[0]; });
        if (uniform && runKinds[0] !== "out") block.classList.add("is-" + runKinds[0]);
        code.insertBefore(block, lines[i]);
        // move lines i..j plus the text nodes between them, and the newline
        // that ends the run, so the block leaves no blank line behind it
        var node = lines[i];
        var end = lines[j];
        while (node) {
          var next = node.nextSibling;
          block.appendChild(node);
          if (node === end) {
            if (next && next.nodeType === 3 && /^\n/.test(next.nodeValue)) {
              if (next.nodeValue === "\n") block.appendChild(next);
              else { next.nodeValue = next.nodeValue.slice(1); block.appendChild(document.createTextNode("\n")); }
            }
            break;
          }
          node = next;
        }
        for (var k = i; k <= j; k++) {
          var line = lines[k];
          line.classList.add("albers-out-line");
          if (kinds[k] !== "out") line.classList.add("is-" + kinds[k]);
          Array.from(line.querySelectorAll(".co")).forEach(function (co) {
            co.classList.add("albers-out");
            if (kinds[k] !== "out") co.classList.add("is-" + kinds[k]);
          });
          wrapPrompt(line);
          stripMark(line);
        }
        i = j + 1;
      }
    });
  }

  /* ---------------------------------------------------------------------
   * Copy buttons (vignettes only; pkgdown ships its own). Copies the source
   * without the "#>" output lines, and announces the result.
   * ------------------------------------------------------------------- */
  var live;
  function announce(msg) {
    if (!live) {
      live = document.createElement("div");
      live.setAttribute("aria-live", "polite");
      live.setAttribute("role", "status");
      live.style.cssText = "position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap;";
      document.body.appendChild(live);
    }
    live.textContent = "";
    setTimeout(function () { live.textContent = msg; }, 30);
  }

  function sourceText(code) {
    var lines = Array.from(code.querySelectorAll(":scope > span:not(.albers-out-block), :scope > .albers-out-block > span"));
    if (!lines.length) return code.innerText;
    return lines
      .filter(function (l) { return !l.classList.contains("albers-out-line"); })
      .map(function (l) { return l.textContent; })
      .join("\n")
      .replace(/\n+$/, "");
  }

  function initCopyButtons() {
    if (isPkgdown()) return;
    if (document.querySelector(".btn-copy-ex, [data-clipboard-copy], div.sourceCode.hasCopyButton")) return;

    document.querySelectorAll("pre.sourceCode > code, div.sourceCode pre > code").forEach(function (code) {
      var pre = code.parentElement;
      if (!pre || pre.querySelector("button.copy-code")) return;
      if (!code.textContent || !code.textContent.trim()) return;
      var btn = document.createElement("button");
      btn.className = "copy-code";
      btn.type = "button";
      btn.setAttribute("aria-label", "Copy code");
      btn.textContent = "copy";
      btn.addEventListener("click", function () {
        var text = sourceText(code);
        var done = function (ok) {
          btn.textContent = ok ? "copied" : (/Mac|iPhone|iPad/.test(navigator.platform) ? "press ⌘C" : "press Ctrl+C");
          btn.classList.add("is-done");
          announce(ok ? "Code copied to clipboard" : "Copy failed; select the code and copy it manually");
          setTimeout(function () { btn.textContent = "copy"; btn.classList.remove("is-done"); }, 1500);
        };
        var fallback = function () {
          try {
            var ta = document.createElement("textarea");
            ta.value = text;
            ta.setAttribute("readonly", "");
            ta.style.cssText = "position:fixed;top:-1000px;opacity:0;";
            document.body.appendChild(ta);
            ta.select();
            var ok = document.execCommand("copy");
            document.body.removeChild(ta);
            btn.focus();
            done(ok);
          } catch (e) { done(false); }
        };
        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(text).then(function () { done(true); }, fallback);
        } else {
          fallback();
        }
      });
      pre.appendChild(btn);
    });
  }

  /* ---------------------------------------------------------------------
   * Figures and tables: numbered captions, scroll wrappers for wide tables
   * ------------------------------------------------------------------- */
  function numberLabel(el, word, n) {
    if (!el || el.querySelector(".albers-num")) return;
    var txt = (el.textContent || "").trim();
    if (!txt || /^(fig(ure)?\.?|table|tab\.?)\s*\d/i.test(txt)) return;
    var span = document.createElement("span");
    span.className = "albers-num";
    span.textContent = word + " " + n;
    el.insertBefore(document.createTextNode(" "), el.firstChild);
    el.insertBefore(span, el.firstChild);
  }

  function initCaptions() {
    var scope = document.querySelector("main#main") || document.body;
    var figs = Array.from(scope.querySelectorAll("div.figure > p.caption, figure > figcaption"));
    figs.forEach(function (cap, i) { numberLabel(cap, "Figure", i + 1); });
    var tabs = Array.from(scope.querySelectorAll("table > caption"));
    tabs.forEach(function (cap, i) { numberLabel(cap, "Table", i + 1); });
  }

  // Long column names (variable_number_12, Sepal.Length) may break after "_"
  // and "." -- on paper, where a wide table can't scroll, that is what lets it
  // fit the page.
  function breakableHeaders(scope) {
    Array.from(scope.querySelectorAll("th")).forEach(function (th) {
      if (th.dataset.albersBreaks || (th.textContent || "").trim().length <= 10) return;
      th.dataset.albersBreaks = "1";
      th.classList.add("albers-th-long");
      var walker = document.createTreeWalker(th, NodeFilter.SHOW_TEXT);
      var nodes = [], n;
      while ((n = walker.nextNode())) nodes.push(n);
      nodes.forEach(function (node) {
        var parts = node.nodeValue.replace(/([_.])/g, "$1\u0000").split("\u0000").filter(function (x) { return x !== ""; });
        if (parts.length < 2) return;
        var frag = document.createDocumentFragment();
        parts.forEach(function (part, i) {
          if (i) frag.appendChild(document.createElement("wbr"));
          frag.appendChild(document.createTextNode(part));
        });
        node.parentNode.replaceChild(frag, node);
      });
    });
  }

  function initTables() {
    var scope = document.querySelector("main#main") || document.body;
    breakableHeaders(scope);
    Array.from(scope.querySelectorAll("table")).forEach(function (t) {
      var row = t.querySelector("tr");
      if (row && row.children.length >= 10) t.classList.add("albers-table-many");
    });
    Array.from(scope.querySelectorAll("table")).forEach(function (t) {
      if (t.closest(".albers-table-scroll, pre, .theme-lab, .gt_table, .dataTables_wrapper, .reactable")) return;
      if (t.parentElement && /table-responsive/.test(t.parentElement.className)) return;
      var wrap = document.createElement("div");
      wrap.className = "albers-table-scroll";
      var cap = t.querySelector(":scope > caption");
      // named only while it scrolls (a region); a plain wrapper is generic
      wrap.dataset.label = cap ? cap.textContent.replace(/\s+/g, " ").trim() : "Table";
      t.parentNode.insertBefore(wrap, t);
      wrap.appendChild(t);
      // A <caption> is only as wide as its table, so a narrow table squeezes
      // it into a column. Set the caption above the scroll region instead;
      // the table keeps an accessible name via aria-labelledby.
      if (cap) {
        var id = "albers-cap-" + Math.random().toString(36).slice(2, 8);
        var p = document.createElement("p");
        p.className = "albers-table-caption";
        p.id = id;
        p.innerHTML = cap.innerHTML;
        wrap.parentNode.insertBefore(p, wrap);
        t.removeChild(cap);
        t.setAttribute("aria-labelledby", id);
      }
    });
  }

  /* Scrollable regions (wide code, wide tables) must be keyboard-reachable;
     regions that do not scroll must not be tab stops. */
  // (on a long page, before it is shown, this waits for idle time: its reads
  // cost a layout of the whole page, and fade cues can come a moment later)
  var shown = false, longPage = null;
  function isLong() {
    if (longPage === null) longPage = document.querySelectorAll("pre > code.sourceCode > span").length > 600;
    return longPage;
  }
  function syncScrollRegions() {
    if (!shown && isLong()) {
      if (!syncScrollRegions.queued) { syncScrollRegions.queued = true; idle(function () { syncScrollRegions.queued = false; syncScrollRegions(); }); }
      return;
    }
    // Measure at the natural measure first: blocks that overflow it are
    // marked .is-wide (on wide screens CSS lets them borrow the margin column).
    var blocks = Array.from(document.querySelectorAll("pre, .albers-table-scroll, .albers-out-block, .albers-math, math:not([display='block'])"));
    blocks.forEach(function (el) {
      el.classList.remove("is-wide");
      if (el.parentElement && el.parentElement.matches("div.sourceCode")) el.parentElement.classList.remove("is-wide");
    });
    blocks.forEach(function (el) {
      if (el.scrollWidth > el.clientWidth + 4 && !el.classList.contains("albers-out-block") && !el.classList.contains("albers-math") && el.tagName.toLowerCase() !== "math") {
        el.classList.add("is-wide");
        if (el.parentElement && el.parentElement.matches("div.sourceCode")) el.parentElement.classList.add("is-wide");
      }
    });
    blocks.forEach(function (el) {
      var scrolls = el.scrollWidth > el.clientWidth + 4;
      el.classList.toggle("is-overflowing", scrolls);
      if (scrolls && !el.dataset.albersScrollWatch) {
        el.dataset.albersScrollWatch = "1";
        el.addEventListener("scroll", function () {
          el.classList.toggle("is-scrolled-end", el.scrollLeft + el.clientWidth >= el.scrollWidth - 2);
        }, { passive: true });
      }
      if (el.tagName.toLowerCase() === "math") {
        // inline math keeps its own semantics: a tab stop only, no region
        if (scrolls) { el.setAttribute("tabindex", "0"); el.dataset.albersTabbed = "1"; }
        else if (el.dataset.albersTabbed) { el.removeAttribute("tabindex"); delete el.dataset.albersTabbed; }
        return;
      }
      if (scrolls) {
        el.setAttribute("tabindex", "0");
        el.setAttribute("role", "region");
        if (!el.hasAttribute("aria-labelledby") && (!el.hasAttribute("aria-label") || /scrolls horizontally/.test(el.getAttribute("aria-label")))) {
          el.setAttribute("aria-label", el.tagName === "PRE" ? "Code (scrolls horizontally)"
            : el.classList.contains("albers-out-block") ? (el.classList.contains("is-msg") ? "Message" : "Output") + " (scrolls horizontally)"
            : el.classList.contains("albers-math") ? (el.dataset.label || "Equation") + " (scrolls horizontally)"
            : (el.dataset.label || "Table") + " (scrolls horizontally)");
        }
      } else if (el.getAttribute("role") === "region") {
        el.removeAttribute("tabindex");
        el.removeAttribute("role");
        if (/scrolls horizontally/.test(el.getAttribute("aria-label") || "")) el.removeAttribute("aria-label");
      }
    });
  }

  // Text for a node with each <math> replaced by its TeX source (the
  // annotation), for names that would otherwise drop the math.
  function textWithMath(node) {
    var c = node.cloneNode(true);
    Array.from(c.querySelectorAll("math")).forEach(function (m) {
      var tex = m.querySelector('annotation[encoding="application/x-tex"]');
      m.replaceWith(document.createTextNode(m.getAttribute("alttext") || (tex ? tex.textContent : m.textContent)));
    });
    return c.textContent.replace(/\s+/g, " ").trim();
  }

  function initMath() {
    // block equations are regions named "Equation n", not raw TeX
    Array.from(document.querySelectorAll('math[display="block"]')).forEach(function (m, i) {
      if (m.parentElement && m.parentElement.classList.contains("albers-math")) return;
      var wrap = document.createElement("div");
      wrap.className = "albers-math";
      wrap.dataset.label = "Equation " + (i + 1);
      m.parentNode.insertBefore(wrap, m);
      wrap.appendChild(m);
    });
    // headings and contents entries keep their inline math in their names
    Array.from(document.querySelectorAll("h1, h2, h3, h4, #TOC a")).forEach(function (el) {
      if (el.querySelector("math") && !el.hasAttribute("aria-label")) el.setAttribute("aria-label", textWithMath(el));
    });
    syncScrollRegions();
  }

  // A callout's run-in label is a bold phrase that opens its first paragraph
  // ("**Tip.** ..."); a bold word later in the sentence is not a label.
  function initCalloutLabels() {
    Array.from(document.querySelectorAll(".callout > p:first-child > strong:first-child")).forEach(function (b) {
      var n = b.previousSibling;
      while (n && n.nodeType === 3 && !/\S/.test(n.nodeValue)) n = n.previousSibling;
      if (!n) b.classList.add("albers-callout-label");
    });
  }

  // a short inline code chip never breaks (no "callout-" / "warning");
  // a long one may, as the last resort
  function initInlineCode() {
    Array.from(document.querySelectorAll(":not(pre) > code")).forEach(function (c) {
      if ((c.textContent || "").length <= 28) c.classList.add("is-short");
    });
  }

  function initA11y() {
    syncScrollRegions();
    // footnote back-links are a bare arrow; name them
    Array.from(document.querySelectorAll(".footnotes li")).forEach(function (li, i) {
      Array.from(li.querySelectorAll("a.footnote-back")).forEach(function (a) {
        if (!a.hasAttribute("aria-label")) a.setAttribute("aria-label", "Back to reference " + (i + 1));
      });
    });
    // widths settle once web fonts and images have loaded
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(syncScrollRegions);
    window.addEventListener("load", syncScrollRegions);
    var t;
    window.addEventListener("resize", function () { clearTimeout(t); t = setTimeout(syncScrollRegions, 150); });
    // Escape closes the phone contents panel
    document.addEventListener("keydown", function (e) {
      if (e.key !== "Escape") return;
      var toc = document.getElementById("TOC");
      if (toc && toc.classList.contains("is-open")) {
        toc.classList.remove("is-open");
        var b = toc.querySelector(".albers-toc-toggle");
        if (b) { b.setAttribute("aria-expanded", "false"); b.focus(); }
      }
    });
    // Vignettes: a skip link to the first section
    if (isVignetteDom() && !document.querySelector(".albers-skip")) {
      var toc = document.getElementById("TOC");
      var first = toc ? toc.nextElementSibling : document.querySelector("body > .section[id]");
      while (first && /^(SCRIPT|STYLE)$/.test(first.tagName)) first = first.nextElementSibling;
      if (first && !first.id) first.id = "albers-content";
      if (first) {
        if (!first.hasAttribute("tabindex")) first.setAttribute("tabindex", "-1");
        var a = document.createElement("a");
        a.className = "albers-skip";
        a.href = "#" + first.id;
        a.textContent = "Skip to content";
        document.body.insertBefore(a, document.body.firstChild);
      }
    }
  }

  /* Print: paper is always light, and folded content is printed open. */
  function initPrint() {
    var saved = null;
    var savedBs = null;
    var opened = [];
    window.addEventListener("beforeprint", function () {
      saved = root.getAttribute("data-albers-theme");
      savedBs = root.getAttribute("data-bs-theme");
      root.setAttribute("data-albers-theme", "light");
      if (savedBs) root.setAttribute("data-bs-theme", "light");
      opened = Array.from(document.querySelectorAll("details:not([open])"));
      opened.forEach(function (d) { d.setAttribute("open", ""); });
    });
    window.addEventListener("afterprint", function () {
      if (saved) root.setAttribute("data-albers-theme", saved);
      if (savedBs) root.setAttribute("data-bs-theme", savedBs);
      opened.forEach(function (d) { d.removeAttribute("open"); });
      opened = [];
    });
  }

  /* ---------------------------------------------------------------------
   * Sidenotes: each footnote is also set in the right margin beside its
   * reference (shown on wide screens only; the numbered list at the end
   * remains for narrow screens and print).
   * ------------------------------------------------------------------- */
  function initSidenotes() {
    var refs = Array.from(document.querySelectorAll("a.footnote-ref[href^='#']"));
    if (!refs.length) return;
    refs.forEach(function (ref) {
      var id = decodeURIComponent(ref.getAttribute("href").slice(1));
      var li = document.getElementById(id);
      if (!li || ref.nextElementSibling && ref.nextElementSibling.classList.contains("albers-sidenote")) return;
      var note = document.createElement("span");
      note.className = "albers-sidenote";
      note.setAttribute("aria-hidden", "true"); // the end-of-page list stays the accessible copy
      var num = document.createElement("span");
      num.className = "albers-sidenote__num";
      num.textContent = (ref.textContent || "").trim();
      note.appendChild(num);
      var body = li.cloneNode(true);
      Array.from(body.querySelectorAll(".footnote-back")).forEach(function (b) { b.parentNode.removeChild(b); });
      Array.from(body.childNodes).forEach(function (n) {
        if (n.nodeType === 1 && n.tagName === "P") {
          var span = document.createElement("span");
          span.className = "albers-sidenote__p";
          span.innerHTML = n.innerHTML;
          note.appendChild(span);
        } else if (n.nodeType === 3 && n.nodeValue.trim()) {
          note.appendChild(document.createTextNode(n.nodeValue));
        }
      });
      Array.from(note.querySelectorAll("a, button, [tabindex]")).forEach(function (x) { x.setAttribute("tabindex", "-1"); });
      ref.parentNode.insertBefore(note, ref.nextSibling);
      // name the reference, and point assistive tech at the note's text
      ref.setAttribute("aria-label", "Footnote " + num.textContent);
      // the description is the note's text alone, without the list item's
      // back-link (a hidden element may still describe)
      var desc = document.createElement("span");
      desc.hidden = true;
      desc.id = "albers-fn-desc-" + num.textContent.replace(/\W+/g, "") + "-" + Math.random().toString(36).slice(2, 6);
      desc.textContent = textWithMath(body);
      ref.parentNode.insertBefore(desc, note.nextSibling);
      ref.setAttribute("aria-describedby", desc.id);
      // While the margin note is on screen, the reference highlights it
      // instead of jumping to the list at the end of the page.
      ref.addEventListener("click", function (e) {
        if (getComputedStyle(note).display === "none") return;
        e.preventDefault();
        note.classList.remove("is-flash");
        void note.offsetWidth;
        note.classList.add("is-flash");
      });
    });
    document.body.classList.add("albers-has-sidenotes");
    // while the margin notes show, the visually hidden list's links leave the tab order
    var list = document.querySelector("body > .footnotes");
    function syncList() {
      if (!list) return;
      // Visually hidden but still read: it is the notes' aria-describedby
      // target (inert would remove it from the accessibility tree). Only its
      // tab stops leave while the margin notes show.
      var hidden = list.getBoundingClientRect().width <= 2;
      Array.from(list.querySelectorAll("a, button, math, [tabindex]")).forEach(function (el) {
        if (hidden) {
          if (!el.hasAttribute("data-albers-tab")) el.setAttribute("data-albers-tab", el.getAttribute("tabindex") || "");
          el.setAttribute("tabindex", "-1");
        } else if (el.hasAttribute("data-albers-tab")) {
          var t = el.getAttribute("data-albers-tab");
          if (t) el.setAttribute("tabindex", t); else el.removeAttribute("tabindex");
          el.removeAttribute("data-albers-tab");
        }
      });
    }
    syncList();
    window.addEventListener("resize", syncList);
  }

  /* Colophon: the page's own swatches and typefaces, set at the foot. */
  function initColophon() {
    if (document.body.classList.contains("albers-no-colophon") || document.querySelector(".albers-colophon")) return;
    var cs = getComputedStyle(document.body);
    var tones = ["--A900", "--A700", "--A500", "--A300", "--sq-outer"];
    var names = ["A900", "A700", "A500", "A300", "comp"];
    var first = function (v) { return (v || "").split(",")[0].replace(/["']/g, "").trim(); };
    var fam = FAMILY_CLASSES.filter(function (f) { return document.body.classList.contains("palette-" + f); })[0] || "red";
    var pre = PRESET_CLASSES.filter(function (f) { return document.body.classList.contains("preset-" + f); })[0] || "homage";
    var foot = document.createElement("footer");
    foot.className = "albers-colophon";
    var strip = document.createElement("div");
    strip.className = "albers-colophon__swatches";
    strip.setAttribute("aria-hidden", "true");
    tones.forEach(function (t, i) {
      var sw = document.createElement("span");
      sw.className = "albers-colophon__swatch";
      sw.style.background = "var(" + t + ")";
      var lab = document.createElement("span");
      lab.textContent = names[i];
      sw.appendChild(lab);
      strip.appendChild(sw);
    });
    var body = first(cs.getPropertyValue("--font-body"));
    var disp = first(cs.getPropertyValue("--font-display"));
    var mono = first(cs.getPropertyValue("--font-mono"));
    var txt = document.createElement("p");
    var faces = body === disp ? body : body + " and " + disp;
    txt.textContent = "Set in " + faces + ", with " + mono + " for code. " +
      "Colour: " + pre.charAt(0).toUpperCase() + pre.slice(1) + " direction, " + fam + " family. Styled with albersdown.";
    foot.appendChild(strip);
    foot.appendChild(txt);
    var notes = document.querySelector("body > .footnotes");
    if (notes && notes.nextSibling) document.body.insertBefore(foot, notes.nextSibling);
    else document.body.appendChild(foot);
  }

  /* ---------------------------------------------------------------------
   * Reading progress (vignettes)
   * ------------------------------------------------------------------- */
  function initProgress() {
    if (document.querySelector(".albers-progress")) return;
    var bar = document.createElement("div");
    bar.className = "albers-progress";
    bar.setAttribute("aria-hidden", "true");
    document.body.appendChild(bar);
    var ticking = false;
    function update() {
      ticking = false;
      var max = document.documentElement.scrollHeight - window.innerHeight;
      var p = max > 0 ? Math.min(1, Math.max(0, window.scrollY / max)) : 0;
      bar.style.setProperty("--albers-progress", p.toFixed(4));
    }
    window.addEventListener("scroll", function () {
      if (!ticking) { ticking = true; nextFrame(update); }
    }, { passive: true });
    window.addEventListener("resize", update);
    update();
  }

  function initDefaults() {
    if (document.body && !hasAny(document.body, "style-", STYLE_CLASSES)) {
      document.body.classList.add("style-minimal");
    }
  }

  /* ---------------------------------------------------------------------
   * Theme Lab
   * ------------------------------------------------------------------- */
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

  function updateLabSummary(labRoot, state) {
    var target = labRoot.querySelector("[data-albers-lab-summary]");
    if (!target) return;
    target.textContent = [
      "family=" + state.family,
      "preset=" + state.preset,
      "style=" + state.style,
      "width=" + (state.width ? state.width + "ch" : "theme default")
    ].join(" | ");
    // the "Copy Into YAML" block follows the choice
    var yaml = document.querySelector("#albers-lab-yaml code");
    if (yaml) {
      yaml.textContent = "output:\n  albersdown::albers_vignette:\n    family: " + state.family +
        "\n    preset: " + state.preset + "\n    style: " + (state.style || "minimal");
    }
  }

  function applyLabState(labRoot, state) {
    applySingleClass(document.body, "palette-", state.family, FAMILY_CLASSES);
    applySingleClass(document.body, "preset-", state.preset, PRESET_CLASSES);

    classes("style-", STYLE_CLASSES).forEach(function (c) {
      document.body.classList.remove(c);
    });
    if (state.style) {
      document.body.classList.add("style-" + state.style);
    }

    // width is a preview only: empty leaves the direction's own measure
    // on <body>: a direction sets its own measure there, over <html>'s
    if (state.width) document.body.style.setProperty("--content", state.width + "ch");
    else document.body.style.removeProperty("--content");
    var widthInput = labRoot.querySelector("[data-albers-control='width']");
    if (widthInput) widthInput.placeholder = state.preset === "interaction" ? "64" : "66";

    applyTheme();

    setSearchParam("family", state.family);
    setSearchParam("preset", state.preset);
    setSearchParam("style", state.style);
    setSearchParam("width", state.width ? String(state.width) : null);

    updateLabSummary(labRoot, state);
    renderAllCompositions();
  }

  function watchThemeClasses() {
    if (!document.body || typeof MutationObserver === "undefined") return;
    var observer = new MutationObserver(function (muts) {
      for (var i = 0; i < muts.length; i++) {
        if (muts[i].attributeName === "class") {
          applyTheme();
          renderAllCompositions();
          break;
        }
      }
    });
    observer.observe(document.body, { attributes: true, attributeFilter: ["class"] });
    // Mirror pkgdown's light-switch into our attribute.
    new MutationObserver(function () {
      if (hasPkgdownLightswitch()) {
        var bs = root.getAttribute("data-bs-theme");
        var want = bs === "dark" ? "dark" : "light";
        if (root.getAttribute("data-albers-theme") !== want) root.setAttribute("data-albers-theme", want);
      }
    }).observe(root, { attributes: true, attributeFilter: ["data-bs-theme"] });
  }

  function initThemeLab() {
    var labRoot = document.querySelector("[data-albers-lab]");
    if (!labRoot) return;

    var controls = {
      family: labRoot.querySelector("[data-albers-control='family']"),
      preset: labRoot.querySelector("[data-albers-control='preset']"),
      style: labRoot.querySelector("[data-albers-control='style']"),
      width: labRoot.querySelector("[data-albers-control='width']")
    };

    if (!controls.family || !controls.preset || !controls.style || !controls.width) return;

    var params = new URLSearchParams(window.location.search);
    // URL values are used only when they name a real choice
    function pick(key, allowed, fallback) {
      var v = (params.get(key) || "").toLowerCase();
      return allowed.indexOf(v) !== -1 ? v : fallback;
    }
    var width = Number(params.get("width") || controls.width.value);
    var initial = {
      family: pick("family", FAMILY_CLASSES, controls.family.value || "red"),
      preset: pick("preset", ["homage", "interaction"], controls.preset.value || "homage"),
      style: pick("style", STYLE_CLASSES, controls.style.value || "minimal"),
      width: width >= 40 && width <= 120 ? width : null
    };

    controls.family.value = initial.family;
    controls.preset.value = initial.preset;
    controls.style.value = initial.style;
    if (initial.width) controls.width.value = String(initial.width);

    applyLabState(labRoot, initial);

    Object.keys(controls).forEach(function (key) {
      controls[key].addEventListener("change", function () {
        applyLabState(labRoot, {
          family: controls.family.value,
          preset: controls.preset.value,
          style: controls.style.value,
          width: Number(controls.width.value) || null
        });
      });
    });
  }

  function safely(fn) {
    try { fn(); } catch (e) { if (window.console) console.warn("albersdown:", e); }
  }

  document.addEventListener("DOMContentLoaded", function () {
    var vignette = isVignetteDom();
    if (vignette) document.body.classList.add("albers-vignette");
    safely(applySiteDefaults);
    safely(initDefaults);
    safely(applyTheme);
    if (vignette) {
      safely(initTitleBlock);
      safely(initToc);
      safely(initProgress);
      safely(initSidenotes);
    }
    if (!vignette) { safely(initSitePlate); safely(initSiteA11y); }
    safely(initFigureTwins);
    safely(initPhoneTwins);
    safely(initWideBlocks);

    safely(initMath);
    safely(initFigureZoom);
    safely(initOutputLines);
    // after output lines are gathered, so only source lines are segmented
    safely(initLineIndents);
    safely(watchSegments);
    safely(initPlainMessages);
    safely(initCalloutLabels);
    safely(initInlineCode);
    safely(initCopyButtons);
    safely(initAnchors);
    safely(initCaptions);
    safely(initTables);
    safely(function () {
      // pkgdown's "Search for" placeholder truncates in a narrow field
      var q = document.querySelector(".navbar input[type='search']");
      if (q && /^Search for/.test(q.getAttribute("placeholder") || "")) q.setAttribute("placeholder", "Search");
    });
    safely(initThemeLab);
    safely(watchThemeClasses);
    safely(renderAllCompositions);
    safely(initPrint);
    safely(initA11y);
    // Inline DOMContentLoaded hooks (older class snippets) run after this
    // listener, in the same event; re-sync the theme and reveal at the next
    // frame, ahead of other queued work (pkgdown's own ready handler held a
    // long article hidden for most of a second). The timer is the fallback
    // where frames do not run (a background tab); `settled` runs once.
    var settled = false;
    function settle() {
      if (settled) return;
      settled = true;
      safely(applySiteDefaults);
      safely(applyTheme);
      safely(function () { initThemeToggle(vignette); });
      if (vignette) safely(initColophon);
      // the browser scrolled to the URL's section before the page grew
      // (colophon, captions); land on it again
      safely(function () {
        if (!location.hash) return;
        var nav = performance.getEntriesByType ? performance.getEntriesByType("navigation")[0] : null;
        if (nav && nav.type !== "navigate") return;
        var el = null;
        try { el = document.getElementById(decodeURIComponent(location.hash.slice(1))); } catch (e) {}
        if (el && window.scrollY > 0) el.scrollIntoView({ block: "start", behavior: "instant" });
      });
      // shown once the fonts are in, but waiting at most 150 ms: the
      // fallback faces are scaled to the web fonts' widths, so a late font
      // rarely reflows the page (a 400 ms wait doubled a vignette's first
      // paint; showing at once let a deck re-wrap from 3 lines to 2)
      var done = false;
      function show() { if (done) return; done = true; reveal(); shown = true; }
      // (not on a long page: its timers run late behind long tasks, and the
      // wait there cost half a second)
      if (!isLong() && document.fonts && document.fonts.status !== "loaded" && document.fonts.ready) {
        document.fonts.ready.then(function () { nextFrame(show); });
        setTimeout(show, 150);
      } else show();
    }
    if (window.requestAnimationFrame) requestAnimationFrame(settle);
    setTimeout(settle, 100);
  });
})();
