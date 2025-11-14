document.addEventListener('DOMContentLoaded', function () {
  // Add copy buttons for standalone R Markdown vignettes (CRAN/offline).
  // Skip if pkgdown has already injected its own clipboard buttons.
  var hasPkgdownCopy =
    document.querySelector('.btn-copy-ex') ||
    document.querySelector('[data-clipboard-copy]') ||
    document.querySelector('div.sourceCode.hasCopyButton');

  if (!hasPkgdownCopy) {
    document.querySelectorAll('pre > code').forEach(function (code) {
      var pre = code.parentElement;
      if (pre.querySelector('button.copy-code')) return; // avoid duplicate buttons we added
      if (!code.textContent || !code.textContent.trim()) return; // skip empty blocks
      var btn = document.createElement('button');
      btn.className = 'copy-code';
      btn.type = 'button';
      btn.setAttribute('aria-label', 'Copy code');
      btn.textContent = 'Copy';
      btn.addEventListener('click', async function () {
        try {
          await navigator.clipboard.writeText(code.innerText);
          btn.textContent = 'Copied';
        } catch (e) {
          // Fallback selection/copy for older browsers
          try {
            var r = document.createRange();
            r.selectNodeContents(code);
            var sel = window.getSelection();
            sel.removeAllRanges();
            sel.addRange(r);
            document.execCommand('copy');
            sel.removeAllRanges();
            btn.textContent = 'Copied';
          } catch (e2) {
            btn.textContent = 'Oops';
          }
        } finally {
          setTimeout(function () { btn.textContent = 'Copy'; }, 1400);
        }
      });
      pre.appendChild(btn);
    });
  }

  // Add visible anchors to h2/h3 (avoid duplicates)
  Array.from(document.querySelectorAll('h2, h3')).forEach(function (h) {
    if (!h.id) return;
    if (h.querySelector('a.anchor')) return;
    var a = document.createElement('a');
    a.href = '#' + h.id;
    a.className = 'anchor';
    a.textContent = '▣';
    a.setAttribute('aria-label', 'Link to this section');
    a.setAttribute('title', 'Link to this section');
    h.appendChild(a);
  });
});
