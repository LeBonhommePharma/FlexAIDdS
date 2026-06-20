// Le Bonhomme Pharma — FlexAID∆S Homepage JS
// Tabs · Copy · Theme · Counter · Drug of day · Mol* hero · Mobile menu

(function () {
  'use strict';

  // ── Tab switching (Usage section) ──────────────────────────
  document.querySelectorAll('.usage-tabs .tab-btn').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var target = btn.getAttribute('aria-controls');
      document.querySelectorAll('.usage-tabs .tab-btn').forEach(function (b) {
        b.classList.remove('active');
        b.setAttribute('aria-selected', 'false');
      });
      document.querySelectorAll('.tab-panel').forEach(function (p) { p.classList.add('hidden'); });
      btn.classList.add('active');
      btn.setAttribute('aria-selected', 'true');
      var panel = document.getElementById(target);
      if (panel) panel.classList.remove('hidden');
    });
  });

  // ── Copy buttons ───────────────────────────────────────────
  document.querySelectorAll('.copy-btn[data-copy]').forEach(function (btn) {
    btn.addEventListener('click', function () {
      navigator.clipboard.writeText(btn.dataset.copy).then(function () {
        var orig = btn.innerHTML;
        btn.innerHTML = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#22D3EE" stroke-width="2.5"><polyline points="20 6 9 17 4 12"/></svg>';
        setTimeout(function () { btn.innerHTML = orig; }, 1500);
      }).catch(function () {});
    });
  });

  // ── Theme toggle ───────────────────────────────────────────
  // Handled by theme.js (anti-flash, localStorage, aria sync).

  // ── Repo stats — seed counters from shared markers ─────────
  function readStat(id, fallback) {
    var el = document.getElementById(id);
    if (!el || !el.textContent) return fallback;
    var n = parseInt(el.textContent.trim(), 10);
    return isNaN(n) ? fallback : n;
  }

  var commitTotal = readStat('stat-commits', 0);
  var langTotal = readStat('stat-langs', 0);

  document.querySelectorAll('[data-count]').forEach(function (el) {
    if (commitTotal > 0) el.dataset.count = String(commitTotal);
  });
  var langDisplay = document.getElementById('stat-langs-display');
  if (langDisplay && langTotal > 0) langDisplay.textContent = String(langTotal);

  function animateCount(el) {
    var target = parseInt(el.dataset.count, 10);
    if (isNaN(target)) return;
    var dur = 1400;
    var start = performance.now();
    function tick(now) {
      var p = Math.min(1, (now - start) / dur);
      var eased = 1 - Math.pow(1 - p, 3);
      el.textContent = Math.round(eased * target);
      if (p < 1) requestAnimationFrame(tick);
      else el.textContent = target;
    }
    requestAnimationFrame(tick);
  }

  if ('IntersectionObserver' in window) {
    var countObs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting) { animateCount(e.target); countObs.unobserve(e.target); }
      });
    }, { threshold: 0.3 });
    document.querySelectorAll('[data-count]').forEach(function (el) { countObs.observe(el); });
  } else {
    document.querySelectorAll('[data-count]').forEach(animateCount);
  }

  // ── Drug of the day (rotation by UTC date) ─────────────────
  // One object drives both the label and the Mol* structure.
  var drugComplexes = [
    { pdb: '1hsg', drug: 'Indinavir',     target: 'HIV-1 protease' },
    { pdb: '3ert', drug: 'Tamoxifen',     target: 'estrogen receptor alpha' },
    { pdb: '1iep', drug: 'Imatinib',      target: 'Abl kinase' },
    { pdb: '1m17', drug: 'Erlotinib',     target: 'EGFR kinase' },
    { pdb: '3nss', drug: 'Oseltamivir',   target: 'influenza neuraminidase' },
    { pdb: '6lu7', drug: 'N3 inhibitor',  target: 'SARS-CoV-2 main protease' },
    { pdb: '4cox', drug: 'Celecoxib',     target: 'cyclooxygenase-2' },
    { pdb: '1hwi', drug: 'Donepezil',     target: 'acetylcholinesterase' },
    { pdb: '2rh1', drug: 'Carazolol',     target: 'beta-2 adrenergic receptor' },
    { pdb: '3htb', drug: 'Dabigatran',    target: 'thrombin' },
    { pdb: '2src', drug: 'Dasatinib',     target: 'Src/Abl kinase' },
    { pdb: '3eml', drug: 'Crizotinib',    target: 'ALK kinase' },
    { pdb: '4dkl', drug: 'Sorafenib',     target: 'RAF kinase' },
    { pdb: '2pgh', drug: 'Flurbiprofen',  target: 'cyclooxygenase' },
    { pdb: '1cbs', drug: 'Retinoic acid', target: 'cellular retinoic acid-binding protein' },
  ];
  var dayIdx = Math.floor(Date.now() / 86400000) % drugComplexes.length;
  var todaysComplex = drugComplexes[dayIdx];

  var label = document.getElementById('drug-of-day-label');
  if (label) {
    label.innerHTML =
      '<strong>' + todaysComplex.drug + '</strong> · ' +
      '<span class="drug-target">' + todaysComplex.target + '</span> ' +
      '<span class="drug-pdb">PDB ' + todaysComplex.pdb.toUpperCase() + '</span>';
  }

  // ── Mobile menu toggle ─────────────────────────────────────
  var mobileToggle = document.querySelector('.mobile-menu-toggle');
  var mainNav = document.querySelector('.main-nav');
  if (mobileToggle && mainNav) {
    mobileToggle.addEventListener('click', function () {
      var isOpen = mobileToggle.getAttribute('aria-expanded') === 'true';
      mobileToggle.setAttribute('aria-expanded', String(!isOpen));
      if (!isOpen) {
        mainNav.style.cssText =
          'display:flex;flex-direction:column;gap:4px;position:absolute;top:56px;left:0;right:0;' +
          'background:rgba(10,14,20,0.98);padding:1rem 1.5rem;border-bottom:1px solid rgba(34,211,238,0.12);z-index:99;';
      } else {
        mainNav.removeAttribute('style');
      }
    });
    mainNav.querySelectorAll('a').forEach(function (a) {
      a.addEventListener('click', function () {
        mobileToggle.setAttribute('aria-expanded', 'false');
        mainNav.removeAttribute('style');
      });
    });
  }

  // ── Sticky header elevation on scroll ─────────────────────
  var header = document.getElementById('site-header');
  if (header) {
    window.addEventListener('scroll', function () {
      header.style.boxShadow = window.scrollY > 10 ? '0 4px 24px rgba(0,0,0,0.4)' : 'none';
    }, { passive: true });
  }

  // ── Mol* viewer — Drug-of-the-Day complex as hero background ─
  var molstarViewer = null;

  function isLightTheme() {
    return document.documentElement.getAttribute('data-theme') === 'light';
  }

  function setPublicationView(viewer) {
    if (!viewer || !viewer.plugin || !viewer.plugin.canvas3d) return;
    try {
      var light = isLightTheme();
      viewer.plugin.canvas3d.setProps({
        renderer: {
          backgroundColor: light ? 0xf8fafc : 0x0a0e14,
          backgroundAlpha: 0,
          ambientIntensity: 0.78,
          lightIntensity: 0.62,
          highlightStrength: 0.35,
        },
        camera: {
          fog: 0,
          clipFar: false,
          clipNear: 0,
        },
        postprocessing: {
          outline: {
            name: 'on',
            params: { scale: 1.15, threshold: 0.22, color: light ? 0x111827 : 0x000000 },
          },
          occlusion: {
            name: 'on',
            params: {
              samples: 48,
              radius: 6,
              bias: 0.8,
              blurKernelSize: 21,
              resolutionScale: 1,
            },
          },
        },
        trackball: {
          animate: { name: 'spin', params: { speed: 0.32 } },
        },
      });
    } catch (e) { /* schema variance across Mol* builds */ }
  }

  function applyDrugOfDayRepresentations(viewer) {
    try {
      var plugin = viewer.plugin;
      var structures = plugin.managers.structure.hierarchy.current.structures;
      if (!structures || !structures.length) return;

      var struct = structures[0];

      plugin.managers.structure.component.clear(struct).then(function () {
        plugin.managers.structure.component.add(
          { structure: struct },
          { type: { name: 'static', params: 'polymer' } }
        ).then(function (polymerComp) {
          if (!polymerComp) return;
          plugin.managers.structure.representation.addRepresentation(polymerComp, {
            type: 'cartoon',
            typeParams: { alpha: 0.42 },
            color: 'chain-id',
          });
        });

        plugin.managers.structure.component.add(
          { structure: struct },
          { type: { name: 'static', params: 'ligand' } }
        ).then(function (ligandComp) {
          if (!ligandComp) return;
          plugin.managers.structure.representation.addRepresentation(ligandComp, {
            type: 'ball-and-stick',
            typeParams: {
              sizeFactor: 0.48,
              sizeAspectRatio: 0.32,
              adjustCylinderLength: true,
              bondScale: 0.55,
              bondSpacing: 0.52,
              linked: true,
              aromaticBonds: true,
              multipleBonds: 'offset',
              includeHydrogens: true,
            },
            color: 'element-symbol',
            colorParams: {
              carbonColor: { name: 'uniform', params: { value: 0xE8E8E8 } },
            },
          });
        });

        plugin.managers.structure.component.add(
          { structure: struct },
          { type: { name: 'static', params: 'branched' } }
        ).then(function (branchedComp) {
          if (!branchedComp) return;
          plugin.managers.structure.representation.addRepresentation(branchedComp, {
            type: 'ball-and-stick',
            typeParams: {
              sizeFactor: 0.44,
              sizeAspectRatio: 0.32,
              adjustCylinderLength: true,
              bondScale: 0.55,
              bondSpacing: 0.52,
              linked: true,
              aromaticBonds: true,
              multipleBonds: 'offset',
              includeHydrogens: true,
            },
            color: 'element-symbol',
          });
        });
      }).then(function () {
        plugin.managers.camera.reset();
      });
    } catch (e) { /* keep default representation if customization fails */ }
  }

  function loadMolstarFallback(done) {
    if (window.molstar) { done(); return; }
    var script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/molstar@4.5.0/build/viewer/molstar.js';
    script.onload = done;
    script.onerror = function () {
      var viewerEl = document.getElementById('molstar-viewer');
      if (viewerEl) viewerEl.classList.add('molstar-unavailable');
    };
    document.head.appendChild(script);
  }

  function initMolstar() {
    if (!window.molstar || !document.getElementById('molstar-viewer')) return;

    molstar.Viewer.create('molstar-viewer', {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowRemoteState: false,
      layoutShowSequence: false,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      viewportShowExpand: false,
      viewportShowSelectionMode: false,
      viewportShowAnimation: false,
      viewportShowControls: false,
      pdbProvider: 'rcsb',
      emdbProvider: 'pdbe',
      canvas3d: {
        transparentBackground: true,
        renderer: { backgroundAlpha: 0 },
        camera: { fog: 0, clipFar: false },
        postprocessing: { outline: { name: 'on' } },
      },
    }).then(function (viewer) {
      molstarViewer = viewer;
      return viewer.loadPdb(todaysComplex.pdb);
    }).then(function () {
      if (!molstarViewer) return;
      return new Promise(function (resolve) {
        setTimeout(function () {
          applyDrugOfDayRepresentations(molstarViewer);
          setPublicationView(molstarViewer);
          resolve();
        }, 1200);
      });
    }).catch(function () {
      var viewerEl = document.getElementById('molstar-viewer');
      if (viewerEl) viewerEl.classList.add('molstar-unavailable');
    });
  }

  function onThemeChange() {
    if (molstarViewer) setPublicationView(molstarViewer);
  }

  if (window.LBPTheme) {
    var origApply = window.LBPTheme.apply;
    window.LBPTheme.apply = function (next) {
      origApply(next);
      onThemeChange();
    };
  } else {
    new MutationObserver(function () { onThemeChange(); })
      .observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
  }

  if (document.readyState === 'complete') {
    loadMolstarFallback(initMolstar);
  } else {
    window.addEventListener('load', function () { loadMolstarFallback(initMolstar); });
  }
})();