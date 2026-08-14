/* ============================================================
   Le Bonhomme Pharma — Canonical theme controller
   - Sets data-theme BEFORE first paint (no flash) when this
     file is loaded synchronously in <head>.
   - Injects the premium pill toggle into every [data-theme-mount].
   - Persists choice in localStorage("lbp-theme").
   - Keeps all toggles + aria state in sync.
   Brand default is dark.
   ============================================================ */
(function () {
  'use strict';

  var KEY = 'lbp-theme';
  var DEFAULT = 'dark';

  function read() {
    try { return localStorage.getItem(KEY); } catch (e) { return null; }
  }
  function write(v) {
    try { localStorage.setItem(KEY, v); } catch (e) { /* private mode */ }
  }

  // ── 1. Anti-flash: apply stored theme immediately (head, pre-paint)
  var theme = read() || DEFAULT;
  document.documentElement.setAttribute('data-theme', theme);

  // ── Toggle markup (single canonical definition) ──────────────
  var SUN =
    '<svg class="lbp-theme-toggle__icon lbp-icon-sun" viewBox="0 0 24 24" fill="none" ' +
    'stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
    '<circle cx="12" cy="12" r="4.2"/>' +
    '<path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/>' +
    '</svg>';
  var MOON =
    '<svg class="lbp-theme-toggle__icon lbp-icon-moon" viewBox="0 0 24 24" fill="none" ' +
    'stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
    '<path d="M21 12.8A9 9 0 1 1 11.2 3 7 7 0 0 0 21 12.8z"/>' +
    '</svg>';

  function makeToggle() {
    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'lbp-theme-toggle';
    btn.setAttribute('role', 'switch');
    btn.setAttribute('title', 'Toggle light / dark theme');
    btn.innerHTML =
      '<span class="lbp-theme-toggle__track">' +
        '<span class="lbp-theme-toggle__thumb">' + SUN + MOON + '</span>' +
      '</span>';
    btn.addEventListener('click', toggle);
    return btn;
  }

  function syncToggles() {
    var isLight = document.documentElement.getAttribute('data-theme') === 'light';
    var label = isLight ? 'Switch to dark mode' : 'Switch to light mode';
    var toggles = document.querySelectorAll('.lbp-theme-toggle');
    for (var i = 0; i < toggles.length; i++) {
      toggles[i].setAttribute('aria-checked', isLight ? 'true' : 'false');
      toggles[i].setAttribute('aria-label', label);
    }
  }

  function apply(next) {
    document.documentElement.setAttribute('data-theme', next);
    write(next);
    syncToggles();
  }

  function toggle() {
    var cur = document.documentElement.getAttribute('data-theme') || DEFAULT;
    apply(cur === 'dark' ? 'light' : 'dark');
  }

  function mountAll() {
    var mounts = document.querySelectorAll('[data-theme-mount]');
    for (var i = 0; i < mounts.length; i++) {
      if (!mounts[i].querySelector('.lbp-theme-toggle')) {
        mounts[i].appendChild(makeToggle());
      }
    }
    syncToggles();
  }

  // ── 2. Inject toggles + enable transitions after first paint ──
  function init() {
    mountAll();

    // React/Babel pages mount their nav after DOMContentLoaded, so keep
    // watching for [data-theme-mount] nodes added later. setTimeout (not
    // rAF) so it still fires in backgrounded / throttled tabs.
    if (window.MutationObserver) {
      var pending = false;
      var obs = new MutationObserver(function () {
        if (pending) return;
        pending = true;
        setTimeout(function () { pending = false; mountAll(); }, 0);
      });
      obs.observe(document.documentElement, { childList: true, subtree: true });
    }

    // Enable cross-fade after first paint so load itself never animates.
    setTimeout(function () {
      document.documentElement.classList.add('theme-anim');
    }, 60);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

  // Expose for any page that wants to react to theme changes.
  window.LBPTheme = { apply: apply, toggle: toggle, current: function () {
    return document.documentElement.getAttribute('data-theme') || DEFAULT;
  } };
})();
