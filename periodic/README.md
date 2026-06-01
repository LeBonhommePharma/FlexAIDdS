# Pharmacological Periodic Table

Self-contained interactive periodic table with pharmacological annotations, designed for deployment at:

**https://thebonhomme.com/periodic**

## Features

- All 118 elements with accurate positioning (including lanthanides/actinides)
- Rich pharmacological profiles for clinically relevant elements:
  - Lithium (mood stabilizer)
  - Fluorine (metabolic blocker in ~25% of drugs)
  - Platinum (cisplatin family — oncology cornerstone)
  - Iodine (thyroid + theranostics)
  - Radium-223 (Xofigo), Lutetium-177 (Pluvicto), Y-90, etc.
  - Gadolinium contrast agents, Bismuth quadruple therapy, Arsenic trioxide (APL), etc.
- Search across symbols, names, drugs, and indications
- Filters: All / Clinically Used / s/p/d/f-block
- Click any tile for detailed mechanism, key drugs, and clinical notes
- Deep linking: `?element=Pt` or `?element=78`
- Keyboard: ⌘K focuses search, Esc clears / closes panel
- Fully matches thebonhomme.com visual language (teal/terra/gold accents, dark scientific palette, IBM Plex + JetBrains Mono)

## Deployment (thebonhomme.com/periodic)

### Option A — Static subfolder (recommended for thebonhomme.com)

1. Copy the entire `periodic/` directory (or just `periodic/index.html`) into the web root of thebonhomme.com under a `periodic/` folder.
2. Ensure your web server serves `periodic/index.html` for `/periodic` and `/periodic/`.

   Example nginx:
   ```nginx
   location /periodic/ {
       alias /path/to/site/periodic/;
       try_files $uri $uri/ /periodic/index.html;
   }
   ```

3. (Optional) Add a redirect for the bare path:
   ```nginx
   location = /periodic {
       return 301 /periodic/;
   }
   ```

### Option B — GitHub Pages / Vercel / Netlify (subpath)

- Deploy the `periodic/` folder as a project or use a monorepo setup with base path `/periodic`.
- For GitHub Pages project site, the URL will be `username.github.io/repo/periodic/`.

### Option C — Full domain redirect note

`entropy.help` currently redirects to `https://thebonhomme.com`. If you want `/periodic` also reachable via `https://entropy.help/periodic`, configure the redirect at the Cloudflare/registrar level to preserve the path, or add a proxy rule.

## Local Development / Testing

```bash
# From repo root
cd periodic
python3 -m http.server 8787
# Open http://localhost:8787
```

Or use any static server:
```bash
npx serve .
```

## Data Sources & Attribution

- Element properties: standard IUPAC / PubChem values
- Pharmacological annotations: synthesized from FDA labels, ATC, clinical literature, and approved drug databases (as of 2026)
- Not for clinical decision support — educational and research tool only

## License

Apache-2.0 (consistent with FlexAIDdS / Le Bonhomme Pharma projects)

---

**Built for https://thebonhomme.com/periodic**  
Le Bonhomme Pharma · Montréal
