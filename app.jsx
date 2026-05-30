// FlexAID∆S Website UI Kit — root app composition.

const { useState: useAppState, useEffect: useAppEffect } = React;

function App() {
  const [active, setActive] = useAppState("hero");

  // Track which section is in view for nav highlighting
  useAppEffect(() => {
    const ids = ["why", "features", "architecture", "install", "benchmarks"];
    const sections = ids.map(id => document.getElementById(id)).filter(Boolean);
    const obs = new IntersectionObserver(entries => {
      const visible = entries
        .filter(e => e.isIntersecting)
        .sort((a, b) => b.intersectionRatio - a.intersectionRatio)[0];
      if (visible) setActive(visible.target.id);
    }, { threshold: [0.2, 0.5, 0.8] });
    sections.forEach(s => obs.observe(s));
    return () => obs.disconnect();
  }, []);

  const jump = (id) => {
    const el = id === "hero" ? document.body : document.getElementById(id);
    if (el) el.scrollIntoView({ behavior: "smooth", block: "start" });
  };

  return (
    <>
      <Nav active={active} onJump={jump} />
      <EntropyMeter />
      <main style={{ paddingTop: "56px", position: "relative", zIndex: 1 }}>
        <HeroSection />
        <div className="section-divider" />
        <WhySection />
        <div className="section-divider" />
        <FeaturesSection />
        <div className="section-divider" />
        <ArchSection />
        <div className="section-divider" />
        <BindingSection />
        <div className="section-divider" />
        <InstallSection />
        <div className="section-divider" />
        <BenchmarksSection />
        <div className="section-divider" />
        <ModulesSection />
        <div className="section-divider" />
        <RepoStatsSection />
      </main>
      <Footer />
    </>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
