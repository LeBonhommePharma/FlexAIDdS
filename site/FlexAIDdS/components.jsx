// FlexAID∆S Website UI Kit — shared primitives.
// Loaded as text/babel; all components export to window for cross-file use.

const { useState, useEffect, useRef } = React;

// ─── Animated conformational-ensemble brand mark (mounts the vanilla engine) ───
function LogoMark({ size = 132, poses = 6, fan = 58, period = 7.2, well = true, className }) {
  const ref = useRef(null);
  useEffect(() => {
    if (!ref.current || !window.mountFlexLogo) return;
    const h = window.mountFlexLogo(ref.current, { size, poses, fan, period, well });
    return () => { h && h.destroy && h.destroy(); };
  }, [size, poses, fan, period, well]);
  return <span ref={ref} className={className} style={{ display: "inline-flex", lineHeight: 0 }} aria-hidden="true" />;
}

// ─── Brand wordmark with animated ∆S ───
function Wordmark({ size = 14 }) {
  return (
    <span className="word" style={{ fontFamily: "var(--font-mono)", fontWeight: 700, color: "var(--teal)", fontSize: size + "px", letterSpacing: "0.02em" }}>
      FlexAID<span className="kw">∆S</span>
    </span>
  );
}

// ─── Section header — eyebrow + h2 ───
function SectionHeader({ eyebrow, children, gradient }) {
  return (
    <>
      <p className="section-label">// {eyebrow}</p>
      <h2 className="section-h2">{children}</h2>
    </>
  );
}

// ─── Feature card variant ───
function FeatureCard({ tint = "teal", title, items }) {
  return (
    <div className={"feature-card " + tint + "-card"}>
      <h3>{title}</h3>
      <ul>{items.map((it, i) => <li key={i}>{it}</li>)}</ul>
    </div>
  );
}

// ─── Architecture pipeline step ───
function ArchStep({ num, title, desc, color = "teal" }) {
  return (
    <div className="arch-step">
      <div className="arch-step-num">Step {num}</div>
      <div className="arch-step-title" style={{ color: "var(--" + color + ")" }}>{title}</div>
      <div className="arch-step-desc">{desc}</div>
      <span className="arch-arrow">→</span>
    </div>
  );
}

// ─── Side-anchored entropy meter ───
function EntropyMeter() {
  const [pct, setPct] = useState(20);
  useEffect(() => {
    const onScroll = () => {
      const max = document.documentElement.scrollHeight - window.innerHeight;
      const p = Math.min(100, Math.max(0, (window.scrollY / Math.max(1, max)) * 100));
      setPct(p);
    };
    window.addEventListener("scroll", onScroll, { passive: true });
    onScroll();
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  // Map scroll position to entropy state
  const state =
    pct < 25 ? { lbl: "UNBOUND", color: "#A78BFA", entropy: (8.5 - pct * 0.05).toFixed(1), dg: "+0.0" } :
    pct < 70 ? { lbl: "ENCOUNTER", color: "#22D3EE", entropy: (5.4 - (pct - 25) * 0.03).toFixed(1), dg: (-(pct - 25) * 0.1).toFixed(1) } :
               { lbl: "BOUND", color: "#FBBF24", entropy: (3.2 - (pct - 70) * 0.02).toFixed(1), dg: (-7.0 - (pct - 70) * 0.05).toFixed(1) };

  return (
    <div className="entropy-meter">
      <span className="em-vert" style={{ color: "#A78BFA" }}>HIGH ENTROPY</span>
      <div className="em-bar">
        <div className="em-fill" style={{ height: pct + "%" }} />
        <div className="em-marker" style={{ top: pct + "%", background: state.color, boxShadow: "0 0 10px " + state.color }} />
      </div>
      <span className="em-vert" style={{ color: "#FBBF24" }}>LOW ENTROPY</span>
      <div>
        <div className="em-v" style={{ color: state.color }}>{state.entropy}</div>
        <div className="em-u">ΔS bits</div>
        <div className="em-v" style={{ color: "var(--gold)", marginTop: "4px" }}>{state.dg}</div>
        <div className="em-u">kcal/mol</div>
      </div>
      <div className="em-state" style={{ color: state.color, background: "color-mix(in srgb, " + state.color + " 12%, transparent)", border: "1px solid color-mix(in srgb, " + state.color + " 40%, transparent)" }}>{state.lbl}</div>
    </div>
  );
}

// ─── Sticky nav ───
function Nav({ active, onJump }) {
  const [open, setOpen] = useState(false);
  const items = ["Why", "Features", "Architecture", "Install", "Benchmarks"];
  const jump = (id) => {
    setOpen(false);
    onJump(id);
  };

  return (
    <nav className="nav">
      <div className="nav-inner">
        <button type="button" className="nav-brand" onClick={() => jump("hero")} aria-label="FlexAID∆S home">
          <LogoMark size={26} poses={4} fan={54} period={6} well={false} />
          <Wordmark size={14} />
        </button>

        <div className="nav-links">
          {items.map(label => {
            const id = label.toLowerCase();
            return (
              <button type="button" key={id}
                 className={"nav-link" + (active === id ? " active" : "")}
                 onClick={() => jump(id)}>{label}</button>
            );
          })}
          <a href="/entropy-driven/" className="nav-link">Entropy</a>
          <a href="/drug-of-the-day/" className="nav-link">Drugs</a>
          <a href="/" className="nav-link nav-home">Home</a>
          <a className="gh nav-link" href="https://github.com/LeBonhommePharma/FlexAIDdS" target="_blank" rel="noreferrer noopener" aria-label="GitHub">
            <svg width="16" height="16" viewBox="0 0 24 24" fill="currentColor" aria-hidden="true"><path d="M12 .5C5.7.5.5 5.7.5 12c0 5.1 3.3 9.4 7.9 10.9.6.1.8-.3.8-.6v-2c-3.2.7-3.9-1.5-3.9-1.5-.5-1.3-1.3-1.7-1.3-1.7-1.1-.7.1-.7.1-.7 1.2.1 1.8 1.2 1.8 1.2 1 1.8 2.8 1.3 3.5 1 .1-.8.4-1.3.7-1.6-2.6-.3-5.3-1.3-5.3-5.7 0-1.3.5-2.3 1.2-3.1-.1-.3-.5-1.5.1-3.1 0 0 1-.3 3.3 1.2a11.5 11.5 0 0 1 6 0C17.3 4.7 18.3 5 18.3 5c.6 1.6.2 2.8.1 3.1.8.8 1.2 1.8 1.2 3.1 0 4.4-2.7 5.4-5.3 5.7.4.4.8 1.1.8 2.2v3.3c0 .3.2.7.8.6 4.6-1.5 7.9-5.8 7.9-10.9C23.5 5.7 18.3.5 12 .5z"/></svg>
          </a>
          <span data-theme-mount></span>
        </div>

        <button type="button"
                className="nav-mobile-btn"
                aria-expanded={open}
                aria-controls="nav-mobile-menu"
                aria-label={open ? "Close menu" : "Open menu"}
                onClick={() => setOpen(v => !v)}>
          <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" aria-hidden="true">
            {open
              ? <><line x1="6" y1="6" x2="18" y2="18"/><line x1="6" y1="18" x2="18" y2="6"/></>
              : <><line x1="4" y1="7" x2="20" y2="7"/><line x1="4" y1="12" x2="20" y2="12"/><line x1="4" y1="17" x2="20" y2="17"/></>}
          </svg>
        </button>
      </div>

      <div id="nav-mobile-menu" className={"nav-mobile-menu" + (open ? " open" : "")}>
        {items.map(label => {
          const id = label.toLowerCase();
          return (
            <button type="button" key={id} className={active === id ? "active" : ""} onClick={() => jump(id)}>{label}</button>
          );
        })}
        <a href="/entropy-driven/" onClick={() => setOpen(false)}>Entropy-Driven</a>
        <a href="/drug-of-the-day/" onClick={() => setOpen(false)}>Drug of the Day</a>
        <a href="/periodic/" onClick={() => setOpen(false)}>Periodic Table</a>
        <a href="/" onClick={() => setOpen(false)}>Home</a>
        <a href="https://github.com/LeBonhommePharma/FlexAIDdS" target="_blank" rel="noreferrer noopener" onClick={() => setOpen(false)}>GitHub</a>
      </div>
    </nav>
  );
}

// ─── Particle drift canvas behind hero ───
function ParticleCanvas() {
  const ref = useRef(null);
  useEffect(() => {
    const c = ref.current;
    if (!c) return;
    const ctx = c.getContext("2d");
    const palette = ["#22D3EE", "#A78BFA", "#FBBF24"];
    let particles = [];
    let raf = 0;
    const resize = () => {
      c.width = c.offsetWidth;
      c.height = c.offsetHeight;
      const n = Math.min(70, Math.floor(c.width * c.height / 18000));
      particles = Array.from({ length: n }, () => ({
        x: Math.random() * c.width, y: Math.random() * c.height,
        vx: (Math.random() - 0.5) * 0.25, vy: (Math.random() - 0.5) * 0.25,
        r: 0.7 + Math.random() * 1.6,
        col: palette[(Math.random() * 3) | 0],
        phase: Math.random() * Math.PI * 2,
      }));
    };
    const step = () => {
      ctx.clearRect(0, 0, c.width, c.height);
      const t = performance.now() / 1000;
      for (const p of particles) {
        p.x += p.vx; p.y += p.vy;
        if (p.x < 0) p.x = c.width; if (p.x > c.width) p.x = 0;
        if (p.y < 0) p.y = c.height; if (p.y > c.height) p.y = 0;
        const a = 0.18 + 0.12 * Math.sin(t * 0.8 + p.phase);
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
        ctx.fillStyle = p.col;
        ctx.globalAlpha = a;
        ctx.fill();
      }
      ctx.globalAlpha = 1;
      raf = requestAnimationFrame(step);
    };
    resize();
    window.addEventListener("resize", resize);
    raf = requestAnimationFrame(step);
    return () => { cancelAnimationFrame(raf); window.removeEventListener("resize", resize); };
  }, []);
  return <canvas ref={ref} className="particle-canvas" aria-hidden="true" />;
}

// ─── Stat that counts up from 0 to its final value on first mount ───
function CountStat({ to, suffix = "", decimals = 0, color, label }) {
  const [v, setV] = useState(0);
  useEffect(() => {
    const dur = 1200;
    const start = performance.now();
    let raf = 0;
    const step = (now) => {
      const p = Math.min(1, (now - start) / dur);
      // expo-out easing
      const eased = 1 - Math.pow(1 - p, 3);
      setV(to * eased);
      if (p < 1) raf = requestAnimationFrame(step);
    };
    raf = requestAnimationFrame(step);
    return () => cancelAnimationFrame(raf);
  }, [to]);
  return (
    <div className="hero-stat">
      <div className="hero-stat-value" style={{ color }}>
        {v.toFixed(decimals)}{suffix}
      </div>
      <div className="hero-stat-label">{label}</div>
    </div>
  );
}

function DrugOfDayBadge() {
  const [drug, setDrug] = React.useState(null);

  React.useEffect(() => {
    function computeTodays() {
      const list = [
        { pdb: '4xp4', drug: 'Cocaine HCl', target: 'dopamine transporter (DAT)', series: '001', href: '/drug-of-the-day/cocaine/' },
        { pdb: '6dzv', drug: 'MDMA',        target: 'serotonin transporter (SERT)', series: '002', href: '/drug-of-the-day/mdma/' },
        { pdb: '6wha', drug: 'DMT',         target: '5-HT\u2082A receptor', series: '003', href: '/drug-of-the-day/dmt/' },
        { pdb: '7wc7', drug: 'Psilocin',    target: '5-HT\u2082A receptor', series: '004', href: '/drug-of-the-day/psilocin/' },
        { pdb: '6wgt', drug: 'LSD',         target: '5-HT\u2082B receptor', series: '005', href: '/drug-of-the-day/lsd/' },
        { pdb: '7xna', drug: 'Amphetamine', target: 'dopamine transporter (DAT)', series: '006', href: '/drug-of-the-day/amphetamine/' },
        { pdb: '8ef5', drug: 'Fentanyl',    target: '\u03bc-opioid receptor', series: '007', href: '/drug-of-the-day/fentanyl/' },
        { pdb: '4djh', drug: 'Salvinorin A', target: '\u03ba-opioid receptor', series: '008', href: '/drug-of-the-day/salvinorin-a/' },
      ];
      const idx = Math.floor(Date.now() / 86400000) % list.length;
      return list[idx];
    }
    const d = (window && window.__TODAYS_DRUG) || computeTodays();
    setDrug(d);
  }, []);

  if (!drug) return null;

  return (
    <a href={drug.href || '/drug-of-the-day/'} className="drug-of-day" aria-label="Drug of the Day">
      <span className="drug-series">#{drug.series}</span>{' '}
      <strong>{drug.drug}</strong> · <span className="drug-target">{drug.target}</span>{' '}
      <span className="drug-pdb">PDB {drug.pdb.toUpperCase()}</span>
    </a>
  );
}

Object.assign(window, { Wordmark, LogoMark, SectionHeader, FeatureCard, ArchStep, EntropyMeter, Nav, ParticleCanvas, CountStat, DrugOfDayBadge });
