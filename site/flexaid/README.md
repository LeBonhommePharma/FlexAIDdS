# /flexaid — Standalone FlexAID∆S Landing Page

This directory serves the dedicated standalone FlexAID∆S marketing / product page at the clean root path:

**https://thebonhomme.com/flexaid**

**Note**: This path is configured as an HTTP redirect (301) to the archive URL via Cloudflare. See the root `CLOUDFLARE_SETUP.md` for the one-time setup the domain owner performs. The HTML file here acts as a fallback.

## Source

The content is the exact standalone HTML file originally located at:

`/Users/lp.more/Downloads/flexaids.html`

It is a self-contained, production-ready landing page featuring:
- Hero with live Mol* 3D molecular viewer background
- Feature grid
- Benchmark results
- Architecture pipeline
- Full dark/light theme support
- Rich structured data (JSON-LD)

## Deployment

Because this lives under `site/flexaid/`, it is automatically included when the "Deploy Site" GitHub Actions workflow runs (`.github/workflows/update-site.yml`).

After merge to `master`, the page will be live at the URL above.

## Local Testing

```bash
# From repo root
python3 -m http.server 8000 --directory site

# Visit: http://localhost:8000/flexaid/
```

## Notes

- The HTML is intentionally kept as a single standalone file (no build step).
- It references external CDNs for Mol* viewer and fonts (same as the original).
- OG metadata currently points to `flexaids.thebonhomme.com` — update if a different canonical domain is chosen.

---

Le Bonhomme Pharma · 2026
