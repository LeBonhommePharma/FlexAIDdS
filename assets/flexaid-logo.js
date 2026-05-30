/* ============================================================
   FlexAID∆S — Dynamic Brand Mark
   Le Bonhomme Pharma · Montréal

   The mark is a living conformational ensemble. A 3-atom ligand
   samples its torsional space as a fan of terra "ghost" poses
   (the −TΔS entropy term). On a slow loop the fan COLLAPSES into
   a single teal bound pose — the minimum on the binding
   free-energy surface (ΔH) — and the bound atoms flare gold at
   the instant of binding (ΔG). Then entropy is released and the
   fan breathes open again.

   ΔG = ΔH − TΔS, rendered as motion.

   Usage:
     mountFlexLogo(elementOrSelector, { size, poses, fan, period, well });
   Returns a handle: { destroy(), el }.
   ============================================================ */
(function () {
  "use strict";

  var NS = "http://www.w3.org/2000/svg";
  var TEAL = "#22D3EE";
  var TERRA = "#A78BFA";
  var GOLD = "#FBBF24";

  // pivot + ligand geometry inside a 100×100 viewBox
  var PX = 50, PY = 52;                 // shared rotation pivot
  var APEX = [50, 24];                  // ligand apex (top atom)
  var L = [28, 64], R = [72, 64];       // two arm atoms

  function lerp(a, b, t) { return a + (b - a) * t; }
  function clamp(v, lo, hi) { return v < lo ? lo : v > hi ? hi : v; }
  function el(tag, attrs) {
    var n = document.createElementNS(NS, tag);
    if (attrs) for (var k in attrs) n.setAttribute(k, attrs[k]);
    return n;
  }

  function buildPose(strokeW, atomR) {
    var g = el("g", {
      "stroke-linecap": "round",
      "stroke-linejoin": "round",
      "fill": "none",
      "stroke-width": strokeW,
      "vector-effect": "non-scaling-stroke"
    });
    var bonds = el("path", { d: "M " + L[0] + " " + L[1] + " L " + APEX[0] + " " + APEX[1] + " L " + R[0] + " " + R[1] });
    g.appendChild(bonds);
    [APEX, L, R].forEach(function (p) {
      g.appendChild(el("circle", { cx: p[0], cy: p[1], r: atomR }));
    });
    return g;
  }

  function mountFlexLogo(target, opts) {
    opts = opts || {};
    var size = opts.size || 132;
    var nGhost = opts.poses || 6;        // ghost poses per side excluded center
    var fan = opts.fan != null ? opts.fan : 58;     // max half-spread in degrees
    var period = (opts.period || 7.2) * 1000;       // full breathe loop (ms)
    var showWell = opts.well !== false;
    var node = typeof target === "string" ? document.querySelector(target) : target;
    if (!node) return null;

    var reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

    var svg = el("svg", {
      viewBox: "0 0 100 100",
      width: size, height: size,
      role: "img", "aria-label": "FlexAID∆S"
    });
    svg.style.display = "block";
    svg.style.overflow = "visible";

    var defs = el("defs", {});
    // free-energy well — radial gradient behind the pivot
    var rg = el("radialGradient", { id: "fxWell" + Math.random().toString(36).slice(2, 7) });
    rg.appendChild(el("stop", { offset: "0%", "stop-color": TEAL, "stop-opacity": "0.22" }));
    rg.appendChild(el("stop", { offset: "55%", "stop-color": TEAL, "stop-opacity": "0.05" }));
    rg.appendChild(el("stop", { offset: "100%", "stop-color": TEAL, "stop-opacity": "0" }));
    defs.appendChild(rg);
    svg.appendChild(defs);

    // ---- layers ----
    var well = el("circle", { cx: PX, cy: PY, r: 34, fill: "url(#" + rg.id + ")" });
    if (showWell) svg.appendChild(well);

    // faint entropy arc the fan sweeps through
    var arc = el("path", {
      fill: "none", stroke: TERRA, "stroke-width": "0.6",
      "stroke-linecap": "round", "stroke-dasharray": "0.5 3", opacity: "0.30"
    });
    svg.appendChild(arc);

    // ghost poses (terra) — symmetric pairs around the bound pose
    var ghosts = [];
    var total = nGhost; // count of ghost layers each carrying a spread fraction
    for (var i = 0; i < total; i++) {
      var g = buildPose(2.0, 2.0);
      g.setAttribute("stroke", TERRA);
      var bonds = g.querySelector("path");
      bonds.setAttribute("stroke", TERRA);
      g.querySelectorAll("circle").forEach(function (c) { c.setAttribute("fill", TERRA); c.setAttribute("stroke", "none"); });
      svg.appendChild(g);
      // spread fraction in (0,1], alternating sides
      var rank = Math.floor(i / 2) + 1;
      var side = i % 2 === 0 ? 1 : -1;
      ghosts.push({ g: g, frac: (rank / Math.ceil(total / 2)) * side });
    }

    // bound pose (teal, on top) — the minimum
    var bound = buildPose(2.6, 2.6);
    bound.setAttribute("stroke", TEAL);
    bound.querySelector("path").setAttribute("stroke", TEAL);
    var boundAtoms = bound.querySelectorAll("circle");
    boundAtoms.forEach(function (c) { c.setAttribute("fill", TEAL); c.setAttribute("stroke", "none"); });
    bound.style.transformBox = "fill-box";
    svg.appendChild(bound);

    node.appendChild(svg);

    function setArc(spreadDeg, entropy) {
      // draw an arc at radius ~ from pivot spanning ±spreadDeg through the apex radius
      var rad = Math.hypot(APEX[0] - PX, APEX[1] - PY);
      var steps = 22, d = "";
      for (var s = 0; s <= steps; s++) {
        var a = lerp(-spreadDeg, spreadDeg, s / steps) * Math.PI / 180;
        // apex direction is straight up from pivot
        var base = Math.atan2(APEX[1] - PY, APEX[0] - PX);
        var ang = base + a;
        var x = PX + rad * Math.cos(ang);
        var y = PY + rad * Math.sin(ang);
        d += (s === 0 ? "M " : "L ") + x.toFixed(2) + " " + y.toFixed(2) + " ";
      }
      arc.setAttribute("d", d);
      arc.setAttribute("opacity", (0.10 + entropy * 0.30).toFixed(3));
    }

    var raf = 0, t0 = null, destroyed = false;

    function frame(now) {
      if (destroyed) return;
      if (t0 == null) t0 = now;
      var t = now - t0;

      // breathe: entropy 0 (collapsed) → 1 (max fan). Asymmetric so the
      // collapse is a quick "snap" and the opening is a slow exhale.
      var phase = (t % period) / period;            // 0..1
      var raw = 0.5 - 0.5 * Math.cos(phase * Math.PI * 2);
      // bias toward staying open, with a sharp dip at the bottom
      var entropy = Math.pow(raw, 0.72);
      var collapse = 1 - entropy;                   // 1 at fully bound

      var curFan = fan * entropy;
      setArc(curFan, entropy);

      // wobble well + glow with collapse (binding event)
      if (showWell) {
        well.setAttribute("r", (28 + collapse * 10).toFixed(2));
        well.setAttribute("opacity", (0.5 + collapse * 0.5).toFixed(3));
      }

      // ghost poses sample the fan, jittering more when entropy is high
      for (var i = 0; i < ghosts.length; i++) {
        var gh = ghosts[i];
        var jitter = Math.sin(t / 1000 * 1.6 + i * 1.7) * 5 * entropy;
        var ang = gh.frac * curFan + jitter;
        gh.g.setAttribute("transform", "rotate(" + ang.toFixed(2) + " " + PX + " " + PY + ")");
        var fade = lerp(0.06, 0.6, entropy) * (1 - 0.28 * Math.abs(gh.frac));
        gh.g.setAttribute("opacity", fade.toFixed(3));
      }

      // bound pose: subtle breathing scale, flares gold at the collapse
      var sc = 1 + collapse * 0.06;
      bound.setAttribute("transform", "translate(" + (PX * (1 - sc)).toFixed(3) + " " + (PY * (1 - sc)).toFixed(3) + ") scale(" + sc.toFixed(3) + ")");
      var glow = 2 + collapse * 9;
      svg.style.filter = "drop-shadow(0 0 " + glow.toFixed(1) + "px rgba(34,211,238," + (0.25 + collapse * 0.45).toFixed(3) + "))";
      // atoms shift teal→gold at the binding instant
      var goldMix = clamp((collapse - 0.55) / 0.45, 0, 1);
      var atomCol = mixHex(TEAL, GOLD, goldMix);
      boundAtoms.forEach(function (c) { c.setAttribute("fill", atomCol); });
      bound.querySelector("path").setAttribute("stroke", mixHex(TEAL, GOLD, goldMix * 0.6));

      raf = requestAnimationFrame(frame);
    }

    if (reduce || opts.static) {
      // static, fully-sampled ensemble
      setArc(fan, 1);
      for (var j = 0; j < ghosts.length; j++) {
        var gj = ghosts[j];
        gj.g.setAttribute("transform", "rotate(" + (gj.frac * fan) + " " + PX + " " + PY + ")");
        gj.g.setAttribute("opacity", (0.34 * (1 - 0.35 * Math.abs(gj.frac))).toFixed(3));
      }
      if (showWell) { well.setAttribute("r", "30"); well.setAttribute("opacity", "0.7"); }
    } else {
      raf = requestAnimationFrame(frame);
    }

    return {
      el: svg,
      destroy: function () { destroyed = true; cancelAnimationFrame(raf); if (svg.parentNode) svg.parentNode.removeChild(svg); }
    };
  }

  function hex(n) { n = Math.round(clamp(n, 0, 255)); var s = n.toString(16); return s.length === 1 ? "0" + s : s; }
  function mixHex(a, b, t) {
    var ar = parseInt(a.slice(1, 3), 16), ag = parseInt(a.slice(3, 5), 16), ab = parseInt(a.slice(5, 7), 16);
    var br = parseInt(b.slice(1, 3), 16), bg = parseInt(b.slice(3, 5), 16), bb = parseInt(b.slice(5, 7), 16);
    return "#" + hex(lerp(ar, br, t)) + hex(lerp(ag, bg, t)) + hex(lerp(ab, bb, t));
  }

  window.mountFlexLogo = mountFlexLogo;
})();
