# /flexaid — FlexAID • Entropy Docking (thebonhomme.com/FlexAID)

Serves the clean dedicated landing at:

**https://thebonhomme.com/flexaid** (and /FlexAID)

This is the focused home for **entropy docking**:

- Prominent "Drug of the Day (DotD)" tag
- **Only the current day's entry appears** (stable daily rotation via day-of-year). No list of everything — just today's.
- Live "Entropy Docking • Tempo of the Beat + Microphone" visual: ensemble states pulse rhythmically (~110 bpm), bottom waveform represents signal. Click the mic to activate live mode (boosts amplitude, jitter, slight tempo shift). Click canvas for tempo bump. Keyboard `m` toggles mic mode.
- Direct ties between thermodynamic entropy collapse and rhythmic/beat metaphors + live signal input (shared kernel language with NATURaL).

The page is intentionally lightweight and self-contained so it renders beautifully at the clean root path.

## Deployment & Cloudflare

- Source lives in `site/flexaid/`.
- CI + GitHub Pages (via the user-site sync) publishes it.
- If a Cloudflare redirect rule still forces `/flexaid` → `/FlexAIDdS/`, remove or adjust it (see `CLOUDFLARE_SETUP.md`). The HTML here is the source of truth for the path.

## Local Testing

```bash
python3 -m http.server 8000 --directory site
# Visit http://localhost:8000/flexaid/
```

## Notes

- Daily DotD logic is pure JS (no server). Uses published entries from the canonical queue.
- Beat + mic visual is canvas + RAF, no external deps.
- Links through to the full `/FlexAIDdS/`, `/entropy-driven/`, and the complete Drug of the Day series.

Le Bonhomme Pharma · 2026
