// FlexAID∆S Website — page-level sections (FlexAIDdS project site).
// Identical to the homepage version except RepoStatsSection reads commit count
// and language count from the hidden <span> elements that update_site_stats.py patches.

const { useState: useStateS, useEffect: useEffectS } = React;

// ─── HERO ───
function HeroSection() {
  return (
    <section id="hero" className="hero">
      <ParticleCanvas />
      <div className="hero-content">
        <p className="hero-tagline-top">Le Bonhomme Pharma · Montréal · <span className="kw">Open Science</span></p>

        <div className="hero-mark">
          <LogoMark size={184} poses={8} fan={60} period={7.2} />
        </div>

        <h1 className="hero-title">
          FlexAID<span className="gold">∆S</span>
        </h1>

        <p className="hero-subtitle">
          <span className="kw">Entropy-Driven</span> <span className="kw">Molecular Docking</span> Engine
        </p>
        <p className="hero-tagline">
          <span className="kw">Genetic algorithms</span> meet <span className="kw">statistical mechanics</span> for real-world drug discovery
        </p>

        <div className="equation">
          <p className="equation-text">
            <span style={{ color: "#FBBF24" }}>ΔG</span>
            <span style={{ color: "#8a93a8" }}> = </span>
            <span style={{ color: "#22D3EE" }}>ΔH</span>
            <span style={{ color: "#8a93a8" }}> − </span>
            <span style={{ color: "#A78BFA" }}>TΔS</span>
          </p>
        </div>

        <div className="hero-stats">
          <CountStat to={0.93} decimals={2} color="#22D3EE" label="PEARSON r (ITC-187) · PRELIMINARY" />
          <CountStat to={1.4}  decimals={1} color="#A78BFA" label="RMSE kcal/mol · PRELIMINARY" />
          <CountStat to={92}   decimals={0} suffix="%" color="#FBBF24" label="BINDING MODE · PRELIMINARY" />
        </div>

        <div className="hero-badges">
          <span className="badge">Apache 2.0</span>
          <span className="badge ga">C++26</span>
          <span className="badge terra">Python ≥ 3.9</span>
          <span className="badge gold">Linux · macOS · Windows</span>
        </div>

        <div className="hero-ctas">
          <a className="btn btn-primary" onClick={() => document.getElementById("install")?.scrollIntoView({ behavior: "smooth" })}>
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
              <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/>
              <polyline points="7 10 12 15 17 10"/>
              <line x1="12" y1="15" x2="12" y2="3"/>
            </svg>
            Get Started
          </a>
          <a className="btn btn-outline" href="https://github.com/LeBonhommePharma/FlexAIDdS" target="_blank" rel="noreferrer noopener">View on GitHub</a>
        </div>
      </div>
    </section>
  );
}

// ─── WHY ───
function WhySection() {
  return (
    <section id="why" className="section">
      <div className="container">
        <SectionHeader eyebrow="why flexaid∆s">
          Why <span className="gradient-tg">FlexAID∆S</span>
        </SectionHeader>
        <div className="note-callout">
          Most <span className="kw">docking engines</span> optimize <span className="kw">enthalpy</span> alone. <strong>FlexAID∆S</strong> adds <span className="kw">conformational entropy</span> via a full <span className="kw">statistical mechanics framework</span> — recovering the correct binding mode <strong>92% of the time</strong> when enthalpy-only scoring fails.
        </div>
      </div>
    </section>
  );
}

// ─── FEATURES ───
function FeaturesSection() {
  return (
    <section id="features" className="section alt">
      <div className="container">
        <SectionHeader eyebrow="capabilities">
          Capabilities <span className="gradient-rg">Grid</span>
        </SectionHeader>
        <div className="features-grid">
          <FeatureCard tint="teal" title="Docking Engine" items={[
            "Genetic algorithm with configurable population, crossover, mutation",
            "Voronoi contact function (CF) for shape complementarity",
            "Dead-end elimination (DEE) torsion pruning",
            "Zero-copy batch scoring via VoronoiCFBatch + OpenMP",
          ]} />
          <FeatureCard tint="terra" title="Thermodynamics" items={[
            "Partition function, free energy, entropy",
            "Shannon entropy + torsional vibrational entropy",
            "Torsional ENCoM (tENCoM) backbone flexibility",
            "Van't Hoff decomposition with ΔCp correction",
          ]} />
          <FeatureCard tint="gold" title="Hardware Acceleration" items={[
            "Unified dispatch: CUDA > Metal > AVX-512 > AVX2 > OpenMP",
            "SIMD primitives for geometric computations",
            "Ultra-fast HPC binaries with LTO + -march=native",
            "Automatic cavity detection (SURFNET + Metal GPU)",
          ]} />
        </div>
      </div>
    </section>
  );
}

// ─── ARCHITECTURE ───
function ArchSection() {
  return (
    <section id="architecture" className="section">
      <div className="container">
        <SectionHeader eyebrow="pipeline">
          Architecture <span className="t-gold">Pipeline</span>
        </SectionHeader>
        <div className="arch-pipeline">
          <ArchStep num="1" title="Input"   color="teal"  desc="PDB receptor + MOL2 ligand. Automatic cavity detection via SURFNET." />
          <ArchStep num="2" title="Sample"  color="terra" desc="GA + DEE torsion pruning. Ring conformer + chirality sampling." />
          <ArchStep num="3" title="Score"   color="gold"  desc="Voronoi CF + ΔS_config via partition function. Free energy F = −kT ln Z." />
          <ArchStep num="4" title="Output"  color="teal"  desc="Ranked poses · ΔG · Van't Hoff decomposition · ITC-validated." />
        </div>
        <div className="arch-sub-grid">
          <div className="arch-sub-card">
            <h4>Flexibility</h4>
            <p><span className="kw">Torsions</span> · Ring conformers · Chirality · <span className="kw">tENCoM backbone</span></p>
          </div>
          <div className="arch-sub-card">
            <h4>Binding Modes</h4>
            <p><span className="kw">Clustering</span> + <span className="kw">ΔG, ΔH, −TΔS</span>, Cv, Boltzmann weights</p>
          </div>
        </div>
      </div>
    </section>
  );
}

// ─── BINDING ANIMATION ───
function BindingSection() {
  const [phase, setPhase] = useStateS(0);
  const phases = [
    {
      lbl: "Diffusion",
      color: "#A78BFA",
      bg: "rgba(167,139,250,0.05)",
      border: "rgba(167,139,250,0.3)",
      desc: <>Drug molecules explore <span className="kw">conformational space</span> freely. High <span className="kw">Shannon entropy</span> reflects many accessible microstates.</>,
      ds: "+8.5", dh: "0.0", dg: "+8.5",
    },
    {
      lbl: "Encounter",
      color: "#22D3EE",
      bg: "rgba(34,211,238,0.05)",
      border: "rgba(34,211,238,0.3)",
      desc: <>The ligand <span className="kw">electrostatically encounters</span> the binding pocket. Translational entropy begins to drop as orientation locks.</>,
      ds: "+4.2", dh: "−3.1", dg: "+1.1",
    },
    {
      lbl: "Binding",
      color: "#FBBF24",
      bg: "rgba(251,191,36,0.05)",
      border: "rgba(251,191,36,0.3)",
      desc: <><span className="kw">Configurational entropy collapses</span> as the ligand locks into the bound pose. Enthalpic interactions dominate ΔG.</>,
      ds: "−5.4", dh: "−12.8", dg: "−7.4",
    },
  ];
  const p = phases[phase];

  return (
    <section id="binding" className="section">
      <div className="container" style={{ maxWidth: "960px" }}>
        <SectionHeader eyebrow="the binding process">
          Drug-Receptor <span className="gradient-rg">Binding</span>
        </SectionHeader>
        <p className="binding-blurb">
          Watch <span className="kw">entropy collapse</span> in real-time as a drug molecule navigates from <span className="kw">chaotic diffusion</span> through <span className="kw">electrostatic encounter</span> to the <span className="kw">bound state</span>.
        </p>

        <div className="binding-canvas-wrap">
          <BindingDiagram phase={phase} />
        </div>

        <div className="binding-controls">
          <div className="binding-phase-btns">
            {phases.map((ph, i) => (
              <button key={i}
                      className={"binding-phase-btn" + (phase === i ? " active-" + (i === 0 ? "diff" : i === 1 ? "enc" : "bind") : "")}
                      onClick={() => setPhase(i)}>
                {i + 1}. {ph.lbl}
              </button>
            ))}
          </div>
        </div>

        <div className="binding-desc" style={{ background: p.bg, border: "1px solid " + p.border }}>
          <p style={{ fontSize: "13px", color: "var(--fg-muted)", lineHeight: 1.6 }}>{p.desc}</p>
          <div className="thermo-row">
            <div><span style={{ color: "var(--fg-muted)" }}>ΔS = </span><span style={{ color: "#A78BFA" }}>{p.ds}</span></div>
            <div><span style={{ color: "var(--fg-muted)" }}>ΔH = </span><span style={{ color: "#22D3EE" }}>{p.dh}</span></div>
            <div><span style={{ color: "var(--fg-muted)" }}>ΔG = </span><span style={{ color: "#FBBF24", fontWeight: 700 }}>{p.dg}</span></div>
          </div>
        </div>
      </div>
    </section>
  );
}

// Static SVG diagram per binding phase
function BindingDiagram({ phase }) {
  const positions = phase === 0
    ? [{x:160,y:60,r:14},{x:340,y:90,r:20},{x:540,y:55,r:11},{x:700,y:120,r:25},{x:480,y:170,r:8}]
    : phase === 1
    ? [{x:220,y:130,r:14},{x:330,y:140,r:18},{x:430,y:130,r:20},{x:520,y:145,r:14}]
    : [{x:450,y:140,r:18}];
  return (
    <svg viewBox="0 0 960 240" style={{ display: "block", width: "100%", height: "240px" }}>
      <defs>
        <radialGradient id="pocket-glow" cx="50%" cy="50%">
          <stop offset="0%" stopColor="#22D3EE" stopOpacity="0.18"/>
          <stop offset="100%" stopColor="#22D3EE" stopOpacity="0"/>
        </radialGradient>
      </defs>
      <ellipse cx="450" cy="140" rx="90" ry="40" fill="url(#pocket-glow)" stroke="#22D3EE" strokeOpacity="0.3" strokeWidth="1" strokeDasharray="3 3"/>
      {phase >= 1 && (
        <path d={phase === 1
          ? "M 220 130 Q 290 100 330 140 Q 380 170 430 130 Q 480 120 520 145"
          : "M 220 130 Q 290 100 330 140 Q 380 170 430 130 Q 460 130 450 140"}
              fill="none" stroke="#FBBF24" strokeWidth="1" strokeOpacity="0.45" strokeDasharray="2 4"/>
      )}
      {positions.map((pos, i) => {
        const isLast = i === positions.length - 1;
        const c = phase === 0 ? "#A78BFA" : phase === 1 ? "#22D3EE" : "#FBBF24";
        return (
          <g key={i} opacity={isLast ? 1 : 0.5}>
            <circle cx={pos.x} cy={pos.y} r={pos.r} fill={c} fillOpacity={isLast ? 0.18 : 0.10}/>
            <g transform={`translate(${pos.x}, ${pos.y}) rotate(${(i * 23) % 60 - 30})`}>
              <path d="M -8 5 L 0 -5 L 8 5" stroke={c} strokeWidth="1.6" fill="none" strokeLinecap="round" strokeLinejoin="round"/>
              <circle cx="-8" cy="5" r="2" fill={c}/>
              <circle cx="0" cy="-5" r="2" fill={c}/>
              <circle cx="8" cy="5" r="2" fill={c}/>
            </g>
          </g>
        );
      })}
    </svg>
  );
}

// ─── MODULES ───
function ModulesSection() {
  const mods = [
    ["tENCoM", <><span className="kw">Torsional ENCoM</span> — backbone normal mode flexibility model</>],
    ["NATURaL", "Co-translational & co-transcriptional assembly"],
    ["ShannonThermoStack", <><span className="kw">Shannon entropy</span> thermodynamic scoring stack</>],
    ["LigandRingFlex", <><span className="kw">Non-aromatic ring conformer sampling</span> + sugar pucker</>],
    ["ChiralCenter", "Explicit R/S stereocenter discrimination"],
    ["CavityDetect", <><span className="kw">SURFNET-based automatic binding site detection</span></>],
  ];
  return (
    <section id="modules" className="section">
      <div className="container">
        <SectionHeader eyebrow="key modules">
          Key <span className="t-gold">Modules</span>
        </SectionHeader>
        <div className="modules-grid">
          {mods.map(([name, desc], i) => (
            <div key={i} className="module-chip">
              <div className="module-name">{name}</div>
              <div className="module-desc">{desc}</div>
            </div>
          ))}
        </div>
      </div>
    </section>
  );
}

// ─── REPO STATS ───
// Commit count and language count are read from hidden <span> elements that
// update_site_stats.py patches on each deploy. The language bar percentages
// stay hardcoded until the next manual refresh of this file.
function RepoStatsSection() {
  const [commits, setCommits] = useStateS(153);
  const [langCount, setLangCount] = useStateS(14);

  useEffectS(() => {
    const c = parseInt(document.getElementById("stat-commits")?.textContent || "153", 10);
    const l = parseInt(document.getElementById("stat-langs")?.textContent || "14", 10);
    if (!isNaN(c)) setCommits(c);
    if (!isNaN(l)) setLangCount(l);
  }, []);

  const langs = [
    ["C++",       57.1, "#f34b7d"],
    ["Python",    21.8, "#3572A5"],
    ["Swift",      6.3, "#F05138"],
    ["TypeScript", 2.6, "#3178c6"],
    ["CMake",      2.5, "#ccc"],
    ["HTML",       1.9, "#e34c26"],
    ["Shell",      1.9, "#89e051"],
    ["Obj-C++",    1.9, "#438eff"],
    ["CUDA",       1.1, "#76B900"],
    ["Other",      2.9, "#555"],
  ];

  return (
    <section className="section alt">
      <div className="container">
        <SectionHeader eyebrow="repository">
          Repository <span className="t-gold">Stats</span>
        </SectionHeader>
        <div className="stats-row">
          <div className="stat-item"><div className="stat-value">{commits}</div><div className="stat-label">Commits</div></div>
          <div className="stat-item"><div className="stat-value">C++26</div><div className="stat-label">Standard</div></div>
          <div className="stat-item"><div className="stat-value">Apache 2.0</div><div className="stat-label">License</div></div>
          <div className="stat-item"><div className="stat-value">{langCount}</div><div className="stat-label">Source Languages</div></div>
        </div>
        <div className="lang-bar">
          {langs.map(([n, w, c]) => (
            <div key={n} style={{ width: w + "%", background: c }} title={`${n} ${w}%`} />
          ))}
        </div>
        <div className="lang-legend">
          {langs.map(([n, w, c]) => (
            <div key={n} className="lang-legend-item">
              <span className="lang-legend-dot" style={{ background: c }} />
              {n} {w}%
            </div>
          ))}
        </div>
      </div>
    </section>
  );
}

// ─── FOOTER ───
function Footer() {
  return (
    <footer className="site-footer">
      <div className="footer-inner">
        <div className="footer-top">
          <div className="footer-brand">
            <LogoMark size={26} poses={4} fan={54} period={6.4} well={false} />
            <Wordmark size={14} />
          </div>
          <div className="footer-links">
            <a href="https://github.com/LeBonhommePharma/FlexAIDdS" target="_blank" rel="noreferrer noopener">GitHub</a>
            <a href="https://x.com/BonhommePharma" target="_blank" rel="noreferrer noopener">@BonhommePharma</a>
            <a href="https://opensource.org/licenses/Apache-2.0" target="_blank" rel="noreferrer noopener">Apache 2.0</a>
          </div>
        </div>
        <div className="footer-bottom">
          <div className="footer-cite">Gaudreault &amp; Najmanovich (2015). J. Chem. Inf. Model. 55(7):1323-36</div>
          <div className="footer-cite">© 2024–2026 LP · Le Bonhomme Pharma · Montréal (Little Burgundy) · Québec · Open Source</div>
          <div className="footer-eq">
            <span className="gold">H(X) = −Σ p(x) log p(x)</span> · Make Entropy Great Again
          </div>
        </div>
      </div>
    </footer>
  );
}

// ─── INSTALL ───
function InstallSection() {
  const [tab, setTab] = useStateS("cli");
  return (
    <section id="install" className="section alt">
      <div className="container" style={{ maxWidth: "896px" }}>
        <SectionHeader eyebrow="get started">
          Installation
        </SectionHeader>
        <div className="install-tabs">
          {[["cli", "CLI Build"], ["python", "Python"]].map(([k, l]) => (
            <button key={k} className={"install-tab" + (tab === k ? " active" : "")} onClick={() => setTab(k)}>{l}</button>
          ))}
        </div>

        {tab === "cli" && (
          <>
            <pre className="code-box">
              <div><span className="plain">$</span> <span className="cmd">git clone</span> https://github.com/LeBonhommePharma/FlexAIDdS.git <span className="plain">&amp;&amp; cd FlexAIDdS</span></div>
              <div><span className="plain">$</span> <span className="cmd">mkdir</span> build <span className="plain">&amp;&amp; cd build</span></div>
              <div><span className="plain">$</span> <span className="cmd">cmake</span> .. -DCMAKE_BUILD_TYPE=Release <span className="plain">&amp;&amp; cmake --build . -j $(nproc)</span></div>
            </pre>
            <table className="cmake-table">
              <thead>
                <tr><th>Option</th><th>Default</th><th>Description</th></tr>
              </thead>
              <tbody>
                <tr><td><code>CMAKE_BUILD_TYPE</code></td><td>Release</td><td>Build type (Debug, Release, RelWithDebInfo)</td></tr>
                <tr><td><code>FLEXAIDS_USE_CUDA</code></td><td>OFF</td><td>Enable CUDA GPU acceleration</td></tr>
                <tr><td><code>FLEXAIDS_USE_METAL</code></td><td>OFF</td><td>Enable Metal GPU acceleration (macOS)</td></tr>
                <tr><td><code>FLEXAIDS_USE_AVX512</code></td><td>OFF</td><td>Enable AVX-512 SIMD</td></tr>
                <tr><td><code>FLEXAIDS_USE_OPENMP</code></td><td>ON</td><td>Enable OpenMP parallel scoring</td></tr>
                <tr><td><code>BUILD_PYTHON_BINDINGS</code></td><td>OFF</td><td>Build Python bindings</td></tr>
              </tbody>
            </table>
          </>
        )}

        {tab === "python" && (
          <pre className="code-box">
            <div><span className="plain">$</span> <span className="cmd">cd</span> python <span className="plain">&amp;&amp; pip install -e .</span></div>
            <div style={{ marginTop: "8px" }}><span className="kw">import</span> flexaidds <span className="kw">as</span> fd</div>
            <div>results = fd.<span className="cmd">dock</span>(receptor=<span className="str">'receptor.pdb'</span>, ligand=<span className="str">'ligand.mol2'</span>, compute_entropy=<span className="kw">True</span>)</div>
          </pre>
        )}
      </div>
    </section>
  );
}

// ─── BENCHMARKS ───
function BenchmarksSection() {
  return (
    <section id="benchmarks" className="section">
      <div className="container">
        <SectionHeader eyebrow="reproducible benchmarks">
          Benchmark <span className="t-gold">Results</span>
        </SectionHeader>
        <p className="binding-blurb">
          Every number generated by <span className="kw">automated benchmark</span> on each commit. Same commit, same container, same seed → <span className="kw">bit-for-bit identical results</span>.
        </p>
        <div className="bench-table-wrap">
          <table className="bench-table">
            <thead>
              <tr>
                <th>Benchmark</th>
                <th className="h">FlexAID∆S</th>
                <th>Vina</th>
                <th>Glide</th>
                <th>rDock</th>
              </tr>
            </thead>
            <tbody>
              <tr><td>ITC-187 ΔG r</td><td className="h">0.93</td><td>0.64</td><td>0.69</td><td>0.61</td></tr>
              <tr><td>ITC-187 RMSE</td><td className="h">1.4</td><td>3.1</td><td>2.9</td><td>3.3</td></tr>
              <tr><td>CASF-2016 Scoring</td><td className="h">0.88</td><td>0.73</td><td>0.78</td><td>0.71</td></tr>
              <tr><td>CNS binding mode rescue</td><td className="h">92%</td><td>58%</td><td>64%</td><td>55%</td></tr>
            </tbody>
          </table>
        </div>
        <p style={{ fontFamily: "var(--font-mono)", fontSize: "10px", color: "var(--fg-muted)", marginTop: "1rem", letterSpacing: "0.08em" }}>
          * Numbers are preliminary; final validation and bootstrap CIs in progress.
        </p>
      </div>
    </section>
  );
}

Object.assign(window, { HeroSection, WhySection, FeaturesSection, ArchSection, BindingSection, BindingDiagram, InstallSection, BenchmarksSection, ModulesSection, RepoStatsSection, Footer });
