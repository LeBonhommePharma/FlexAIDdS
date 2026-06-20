# /periodic — Pharmacological Periodic Table

This directory is served at the clean root path **https://thebonhomme.com/periodic** (and via any redirects such as entropy.help).

**Current requirement**: Cloudflare (free) URL rewrites are used in front of GitHub Pages to remove the `/FlexAIDdS/` project-site prefix. See the committed root file `CLOUDFLARE_SETUP.md` for the exact one-time steps the domain owner must perform.

## How it is published

The existing GitHub Actions workflow (`.github/workflows/update-site.yml`) deploys the entire `site/` directory as the GitHub Pages root. Any change under `site/**` (including this folder) triggers a new deployment.

Once the Cloudflare rewrites are active, content appears at the clean URL with no further workflow changes required.

## Contents

- `index.html` — Self-contained, Tailwind-CDN-powered interactive periodic table with deep pharmacological annotations for clinically relevant elements (Li, Pt, F, I, Ra-223, Lu-177, Gd agents, Bi, As₂O₃, etc.).
- Fully matches thebonhomme.com / FlexAID∆S design system.

## Local testing

```bash
# From repo root
python3 -m http.server 8000 --directory site
# Then visit http://localhost:8000/periodic/
```

## Adding / updating content

Edit `site/periodic/index.html` directly. The file is deliberately standalone (no build step required) so it can be iterated quickly while still living inside the main site deployment pipeline.

---

Le Bonhomme Pharma · 2026
