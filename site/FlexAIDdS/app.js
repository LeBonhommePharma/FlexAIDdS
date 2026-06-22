// Le Bonhomme Pharma — FlexAID∆S Homepage JS
// Tabs · Copy · Theme · Counter · Drug of day · Mol* hero · Mobile menu

(function () {
  'use strict';

  var MOLSTAR_BUILD = '20260620-theme';

  function installMolstarNetworkGuards() {
    if (window.__FLEXAID_MOLSTAR_GUARDS__) return;
    window.__FLEXAID_MOLSTAR_GUARDS__ = true;
    var origFetch = window.fetch;
    if (!origFetch) return;
    window.fetch = function (input, init) {
      var url = typeof input === 'string' ? input : (input && input.url) || '';
      if (url.indexOf('molstarvolseg.ncbr.muni.cz') !== -1) {
        return Promise.resolve(new Response('{}', {
          status: 200,
          headers: { 'Content-Type': 'application/json' },
        }));
      }
      return origFetch.apply(this, arguments);
    };
  }

  function bindMolstarRejectionGuard() {
    if (window.__FLEXAID_MOLSTAR_REJECTION_GUARD__) return;
    window.__FLEXAID_MOLSTAR_REJECTION_GUARD__ = true;
    window.addEventListener('unhandledrejection', function (e) {
      var reason = e.reason;
      var msg = reason && (reason.message || String(reason)) || '';
      if (/molstarvolseg|multiScale|is not iterable|Load failed/i.test(msg)) {
        e.preventDefault();
        molstarLog('warn', 'suppressed non-fatal Mol* rejection', reason);
      }
    });
  }

  installMolstarNetworkGuards();

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

  // ── Drug of the Day (published #001–#003; queue #004–#008 in drug-of-the-day/queue.json) ──
  var drugComplexes = [
    { pdb: '4xp4', drug: 'Cocaine HCl', target: 'dopamine transporter (DAT)', series: '001', href: '/drug-of-the-day/cocaine/' },
    { pdb: '6dzv', drug: 'MDMA',        target: 'serotonin transporter (SERT)', series: '002', href: '/drug-of-the-day/mdma/' },
    { pdb: '6wha', drug: 'DMT',         target: '5-HT\u2082A receptor', series: '003', href: '/drug-of-the-day/dmt/' },
    { pdb: '7wc7', drug: 'Psilocin',    target: '5-HT\u2082A receptor', series: '004', href: '/drug-of-the-day/psilocin/' },
    { pdb: '6wgt', drug: 'LSD',         target: '5-HT\u2082B receptor', series: '005', href: '/drug-of-the-day/lsd/' },
    { pdb: '7xna', drug: 'Amphetamine', target: 'dopamine transporter (DAT)', series: '006', href: '/drug-of-the-day/amphetamine/' },
    { pdb: '8ef5', drug: 'Fentanyl',    target: '\u03bc-opioid receptor', series: '007', href: '/drug-of-the-day/fentanyl/' },
    { pdb: '4djh', drug: 'Salvinorin A', target: '\u03ba-opioid receptor', series: '008', href: '/drug-of-the-day/salvinorin-a/' },
  ];
  var dayIdx = Math.floor(Date.now() / 86400000) % drugComplexes.length;
  var todaysComplex = drugComplexes[dayIdx];

  var label = document.getElementById('drug-of-day-label');
  if (label) {
    label.innerHTML =
      '<span class="drug-series">#' + todaysComplex.series + '</span> ' +
      '<strong>' + todaysComplex.drug + '</strong> · ' +
      '<span class="drug-target">' + todaysComplex.target + '</span> ' +
      '<span class="drug-pdb">PDB ' + todaysComplex.pdb.toUpperCase() + '</span>';
  }
  var drugLink = document.querySelector('.drug-of-day');
  if (drugLink && todaysComplex.href) drugLink.setAttribute('href', todaysComplex.href);

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
  var molstarResizeObserver = null;

  function isIOSLike() {
    var ua = navigator.userAgent || '';
    if (/iPad|iPhone|iPod/i.test(ua)) return true;
    return navigator.platform === 'MacIntel' && navigator.maxTouchPoints > 1;
  }

  function molstarLog(level, message, detail) {
    var fn = console[level] || console.log;
    if (detail !== undefined) fn.call(console, '[molstar] ' + message, detail);
    else fn.call(console, '[molstar] ' + message);
  }

  bindMolstarRejectionGuard();

  function isLightTheme() {
    return document.documentElement.getAttribute('data-theme') === 'light';
  }

  function markViewerReady() {
    var el = document.getElementById('molstar-viewer');
    if (!el) return;
    el.classList.remove('molstar-unavailable');
    if (el.querySelector('canvas')) el.classList.add('molstar-ready');
  }

  function markViewerUnavailable(reason, detail) {
    molstarLog('warn', reason, detail);
    var el = document.getElementById('molstar-viewer');
    if (el) el.classList.add('molstar-unavailable');
  }

  function hasRenderableCanvas() {
    var el = document.getElementById('molstar-viewer');
    var canvas = el && el.querySelector('canvas');
    return !!(canvas && canvas.width > 0 && canvas.height > 0);
  }

  function recoverViewerIfCanvasReady(reason, detail) {
    if (hasRenderableCanvas()) {
      molstarLog('warn', reason + ' — canvas rendered, keeping viewer visible', detail);
      markViewerReady();
      if (molstarViewer && molstarViewer.plugin && molstarViewer.plugin.canvas3d) {
        molstarViewer.plugin.canvas3d.requestDraw();
      }
      return true;
    }
    markViewerUnavailable(reason, detail);
    return false;
  }

  function waitForStructures(plugin, maxMs) {
    maxMs = maxMs || 15000;
    var started = Date.now();
    return new Promise(function (resolve) {
      function poll() {
        var structures = plugin.managers.structure.hierarchy.current.structures;
        if (structures && structures.length) {
          resolve(structures);
          return;
        }
        if (Date.now() - started >= maxMs) {
          resolve(null);
          return;
        }
        setTimeout(poll, 250);
      }
      poll();
    });
  }

  function viewerThemePalette() {
    var light = isLightTheme();
    return {
      backgroundColor: light ? 0xffffff : 0x0a0e14,
      outlineColor: light ? 0x000000 : 0xffffff,
      backgroundAlpha: 1,
    };
  }

  function safePostprocessing(outlineOn) {
    var palette = viewerThemePalette();
    var ios = isIOSLike();
    var outline = {
      name: outlineOn ? 'on' : 'off',
      params: outlineOn ? {
        scale: ios ? 1.05 : 1.15,
        threshold: ios ? 0.26 : 0.22,
        color: palette.outlineColor,
      } : {},
    };
    return { outline: outline, occlusion: { name: 'off', params: {} } };
  }

  function publicationViewProps(light, tier) {
    var palette = viewerThemePalette();
    var ios = isIOSLike();
    var props = {
      renderer: {
        backgroundColor: palette.backgroundColor,
        backgroundAlpha: palette.backgroundAlpha,
        ambientIntensity: ios ? 0.85 : 0.78,
        lightIntensity: ios ? 0.68 : 0.62,
        highlightStrength: ios ? 0.28 : 0.35,
      },
      camera: {
        fog: 0,
        clipFar: false,
        clipNear: 0,
      },
      trackball: {
        animate: { name: 'spin', params: { speed: ios ? 0.24 : 0.32 } },
      },
    };

    if (tier === 'outline') {
      props.postprocessing = safePostprocessing(true);
    } else if (tier === 'basic') {
      props.postprocessing = safePostprocessing(false);
    }

    return props;
  }

  function setPublicationView(viewer) {
    if (!viewer || !viewer.plugin || !viewer.plugin.canvas3d) return;
    var light = isLightTheme();
    var canvas3d = viewer.plugin.canvas3d;
    var tiers = ['outline', 'basic'];

    function applyTier(index) {
      if (index >= tiers.length) return;
      try {
        canvas3d.setProps(publicationViewProps(light, tiers[index]));
        canvas3d.requestDraw();
      } catch (e) {
        molstarLog('warn', 'publication view tier failed: ' + tiers[index], e);
        applyTier(index + 1);
      }
    }

    applyTier(0);
  }

  function ensureViewerLayout(container) {
    var hero = container.closest('.hero');
    var height = 480;
    if (hero) {
      var rect = hero.getBoundingClientRect();
      if (rect.height > 0) height = Math.round(rect.height);
    }
    container.style.minHeight = height + 'px';
    container.style.width = '100%';
  }

  function bindViewerResize(viewer, container) {
    function refresh() {
      ensureViewerLayout(container);
      try {
        if (viewer.plugin && viewer.plugin.layout && viewer.plugin.layout.events) {
          viewer.plugin.layout.events.updated.next();
        }
        if (viewer.plugin && viewer.plugin.canvas3d) {
          viewer.plugin.canvas3d.requestResize();
          viewer.plugin.canvas3d.requestDraw();
        }
      } catch (e) { /* Safari layout refresh is best-effort */ }
    }

    if (molstarResizeObserver) molstarResizeObserver.disconnect();
    if ('ResizeObserver' in window) {
      molstarResizeObserver = new ResizeObserver(function () { refresh(); });
      molstarResizeObserver.observe(container);
      if (container.parentElement) molstarResizeObserver.observe(container.parentElement);
    }
    window.addEventListener('resize', refresh, { passive: true });
    window.addEventListener('orientationchange', function () { setTimeout(refresh, 250); }, { passive: true });
    if (window.visualViewport) {
      window.visualViewport.addEventListener('resize', refresh, { passive: true });
    }
    refresh();
  }

  function whenLayoutReady(container, cb) {
    ensureViewerLayout(container);
    requestAnimationFrame(function () {
      requestAnimationFrame(cb);
    });
  }

  function bindRenderErrorRecovery(viewer) {
    try {
      var canvas3d = viewer.plugin.canvas3d;
      if (!canvas3d || canvas3d.__flexaidRenderGuard) return;
      canvas3d.__flexaidRenderGuard = true;
      var origUpdate = canvas3d.update.bind(canvas3d);
      var recovered = false;
      canvas3d.update = function () {
        try {
          return origUpdate();
        } catch (e) {
          var msg = e && (e.message || String(e)) || '';
          if (!recovered && /multiScale/i.test(msg)) {
            recovered = true;
            molstarLog('warn', 'render loop error — disabling postprocessing', e);
            canvas3d.setProps({ postprocessing: safePostprocessing(false) });
            try { return origUpdate(); } catch (e2) { /* swallow repeat */ }
          }
        }
      };
    } catch (e) { /* optional guard */ }
  }

  function bindWebGLContextGuards(viewer) {
    try {
      var canvas = viewer.plugin.canvas3d.canvas;
      if (!canvas) return;
      canvas.addEventListener('webglcontextlost', function (e) {
        e.preventDefault();
        molstarLog('warn', 'WebGL context lost — attempting restore');
      }, false);
      canvas.addEventListener('webglcontextrestored', function () {
        molstarLog('warn', 'WebGL context restored — redrawing');
        setPublicationView(viewer);
        if (viewer.plugin.canvas3d) viewer.plugin.canvas3d.requestDraw();
      }, false);
    } catch (e) { /* optional guard */ }
  }

  function ballAndStickParams() {
    return {
      sizeFactor: 0.48,
      sizeAspectRatio: 0.32,
      adjustCylinderLength: true,
      bondScale: 0.55,
      bondSpacing: 0.52,
      linked: true,
      aromaticBonds: true,
      multipleBonds: 'offset',
      includeHydrogens: true,
    };
  }

  function addRepresentation(plugin, component, spec) {
    if (!component) return Promise.resolve(false);
    return plugin.managers.structure.representation.addRepresentation(component, spec)
      .then(function () { return true; })
      .catch(function (err) {
        molstarLog('warn', 'representation add failed', err);
        return false;
      });
  }

  function applyDrugOfDayRepresentations(viewer) {
    var plugin = viewer.plugin;
    return waitForStructures(plugin).then(function (structures) {
      if (!structures || !structures.length) {
        molstarLog('warn', 'no structures loaded — keeping default representation');
        return false;
      }

      var struct = structures[0];

      function addCustomRepresentations() {
        var polymerOk = false;
        var ligandOk = false;

        return plugin.managers.structure.component.add(
          { structure: struct },
          { type: { name: 'static', params: 'polymer' } }
        ).then(function (polymerComp) {
          return addRepresentation(plugin, polymerComp, {
            type: 'cartoon',
            typeParams: { alpha: isIOSLike() ? 0.5 : 0.42 },
            color: 'chain-id',
          }).then(function (ok) { polymerOk = ok; return ok; });
        }).then(function () {
          return plugin.managers.structure.component.add(
            { structure: struct },
            { type: { name: 'static', params: 'ligand' } }
          );
        }).then(function (ligandComp) {
          return addRepresentation(plugin, ligandComp, {
            type: 'ball-and-stick',
            typeParams: ballAndStickParams(),
            color: 'element-symbol',
            colorParams: {
              carbonColor: { name: 'uniform', params: { value: 0xE8E8E8 } },
            },
          }).then(function (ok) { ligandOk = ok; return ok; });
        }).then(function () {
          return plugin.managers.structure.component.add(
            { structure: struct },
            { type: { name: 'static', params: 'branched' } }
          );
        }).then(function (branchedComp) {
          return addRepresentation(plugin, branchedComp, {
            type: 'ball-and-stick',
            typeParams: ballAndStickParams(),
            color: 'element-symbol',
          });
        }).then(function () {
          plugin.managers.camera.reset();
          return polymerOk || ligandOk;
        });
      }

      return plugin.managers.structure.component.clear(struct).then(function () {
        return addCustomRepresentations();
      }).catch(function (err) {
        molstarLog('warn', 'clear/default customize failed — adding on top of defaults', err);
        return addCustomRepresentations();
      }).catch(function (err) {
        molstarLog('warn', 'custom representations failed — keeping defaults', err);
        plugin.managers.camera.reset();
        return false;
      });
    });
  }

  function loadStructure(viewer, pdbId) {
    var plugin = viewer.plugin;
    var id = pdbId.toLowerCase();

    return viewer.loadPdb(id).catch(function (err) {
      molstarLog('warn', 'loadPdb failed for ' + id + ', trying RCSB PDB', err);
      return plugin.loadStructureFromUrl(
        'https://files.rcsb.org/download/' + id + '.pdb',
        'pdb',
        { label: id.toUpperCase() }
      );
    }).catch(function (err) {
      molstarLog('warn', 'RCSB PDB failed for ' + id + ', trying PDBe bcif', err);
      return plugin.loadStructureFromUrl(
        'https://www.ebi.ac.uk/pdbe/entry-files/download/' + id + '.bcif',
        'bcif',
        { label: id.toUpperCase() }
      );
    }).then(function () { return true; })
      .catch(function (err) {
        molstarLog('error', 'all structure loaders failed for ' + id, err);
        return false;
      });
  }

  function loadMolstarScript(done) {
    if (window.molstar) {
      done();
      return;
    }
    var script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/molstar@4.5.0/build/viewer/molstar.js';
    script.onload = done;
    script.onerror = function () {
      var fallback = document.createElement('script');
      fallback.src = 'https://unpkg.com/molstar@4.5.0/build/viewer/molstar.js';
      fallback.onload = done;
      fallback.onerror = function () {
        markViewerUnavailable('Mol* script failed to load from CDN');
      };
      document.head.appendChild(fallback);
    };
    document.head.appendChild(script);
  }

  function initMolstar() {
    var container = document.getElementById('molstar-viewer');
    if (!container) return;
    if (!window.molstar) {
      markViewerUnavailable('Mol* library not available');
      return;
    }

    whenLayoutReady(container, function () {
      var ios = isIOSLike();
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
        volumesAndSegmentationsDefaultServer: '',
        volumeStreamingDisabled: true,
        canvas3d: {
          transparentBackground: false,
          renderer: {
            backgroundColor: viewerThemePalette().backgroundColor,
            backgroundAlpha: 1,
            pixelRatio: ios ? Math.min(window.devicePixelRatio || 1, 2) : (window.devicePixelRatio || 1),
          },
          camera: { fog: 0, clipFar: false },
          postprocessing: safePostprocessing(true),
          trackball: {
            noScroll: true,
            noRotate: false,
            noPan: true,
            noZoom: true,
          },
        },
      }).then(function (viewer) {
        molstarViewer = viewer;
        bindViewerResize(viewer, container);
        bindRenderErrorRecovery(viewer);
        bindWebGLContextGuards(viewer);
        return loadStructure(viewer, todaysComplex.pdb);
      }).then(function (loaded) {
        if (!loaded || !molstarViewer) {
          markViewerUnavailable('structure failed to load for PDB ' + todaysComplex.pdb);
          return;
        }
        setPublicationView(molstarViewer);
        var applyReps = (window.MolstarDrugViewer && window.MolstarDrugViewer.applyPublicationRepresentations)
          ? window.MolstarDrugViewer.applyPublicationRepresentations.bind(window.MolstarDrugViewer)
          : applyDrugOfDayRepresentations;
        return applyReps(molstarViewer).catch(function (repErr) {
          molstarLog('warn', 'representation customize failed — keeping default render', repErr);
        });
      }).then(function () {
        if (!molstarViewer) return;
        markViewerReady();
        if (molstarViewer.plugin && molstarViewer.plugin.canvas3d) {
          molstarViewer.plugin.canvas3d.requestDraw();
        }
      }).catch(function (err) {
        recoverViewerIfCanvasReady('viewer initialization failed', err);
      });
    });
  }

  function onThemeChange() {
    if (!molstarViewer) return;
    setPublicationView(molstarViewer);
    if (molstarViewer.plugin && molstarViewer.plugin.canvas3d) {
      molstarViewer.plugin.canvas3d.requestDraw();
    }
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

  function bootMolstar() {
    var attempts = 0;
    function tryBoot() {
      if (window.molstar) {
        initMolstar();
        return;
      }
      if (attempts++ < 40) {
        setTimeout(tryBoot, 100);
        return;
      }
      loadMolstarScript(initMolstar);
    }
    tryBoot();
  }

  function startWhenViewerReady() {
    if (document.getElementById('molstar-viewer')) {
      bootMolstar();
      return;
    }
    var obs = new MutationObserver(function () {
      if (document.getElementById('molstar-viewer')) {
        obs.disconnect();
        bootMolstar();
      }
    });
    var root = document.body || document.documentElement;
    obs.observe(root, { childList: true, subtree: true });
    setTimeout(function () { obs.disconnect(); }, 30000);
  }

  if (document.readyState === 'complete') {
    startWhenViewerReady();
  } else {
    window.addEventListener('load', startWhenViewerReady);
  }

  window.__FLEXAID_MOLSTAR_BUILD__ = MOLSTAR_BUILD;
  window.bootFlexaidMolstar = bootMolstar;
})();