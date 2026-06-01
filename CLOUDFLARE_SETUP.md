# Cloudflare Setup for Clean Root Paths (thebonhomme.com/periodic, /flexaid, /entropy-driven)

## Why this is needed
The "Deploy Site" workflow (`.github/workflows/update-site.yml`) uploads the entire `site/` folder as the GitHub Pages artifact.

Because this is a **project site** (not a user/org `*.github.io` root site), GitHub serves everything under the repo slug prefix:
- Deployed content lives at `https://thebonhomme.com/FlexAIDdS/periodic` etc.
- Confirmed in raw workflow logs: "Evaluated environment url: http://thebonhomme.com/FlexAIDdS/"

You requested clean apex paths:
- https://thebonhomme.com/periodic (full interactive pharmacological periodic table)
- https://thebonhomme.com/flexaid (redirect to the archive pre-redesign page)
- https://thebonhomme.com/entropy-driven (standalone page)

**Solution chosen**: Put Cloudflare (free tier) in front of the custom domain and use **URL Rewrite Transform Rules** (plus one Redirect Rule). This requires **zero changes** to the existing deployment pipeline, `site/` structure, or future PRs.

Repo-side preparation (already completed and pushed to master):
- Hygiene improvements to `.gitignore` (prevents future untracked campaign/monitor script merge blocks).
- Updated the limitation note in `site/periodic/index.html` (and the standalone enhanced copy) to accurately describe the GitHub Pages + Cloudflare setup.

## Your required actions (registrar + Cloudflare dashboard only)

### 1. Add the domain to Cloudflare
1. Go to https://dash.cloudflare.com
2. "Add a Site" → enter `thebonhomme.com`
3. Choose the **Free** plan.
4. Cloudflare will scan your existing DNS records.

### 2. Change nameservers at your domain registrar
This step is **mandatory** for Cloudflare proxying (orange cloud) and Transform Rules to work on an apex domain.

- Cloudflare will display two nameservers (e.g. `ada.ns.cloudflare.com` and another).
- Log in to the registrar where you purchased `thebonhomme.com`.
- Replace the current nameservers with Cloudflare's two nameservers.
- Save. Propagation usually takes 5–60 minutes.

**Do not** use "CNAME setup" or partial setups if you want full control (Transform Rules + apex support). Full nameserver delegation is the standard and most reliable path.

### 3. Configure DNS records in Cloudflare
After nameserver propagation, go to the DNS tab for the zone.

**Apex domain (`@` / thebonhomme.com)** — create four A records, **all Proxied** (orange cloud icon):
- `185.199.108.153`
- `185.199.109.153`
- `185.199.110.153`
- `185.199.111.153`

**www (recommended)**:
- Type: CNAME
- Name: `www`
- Target: `LeBonhommePharma.github.io` (confirm the exact value shown in your GitHub repo's Pages settings)
- Proxy status: Proxied (orange)

Delete or disable any old conflicting A/CNAME records that pointed directly elsewhere.

### 4. Confirm / set the custom domain in GitHub Pages
1. Go to the repo Settings → Pages.
2. Under "Custom domain", ensure it says `thebonhomme.com`.
3. Check "Enforce HTTPS" if the option appears (it may take a few minutes after DNS is correct).

The domain verification should already be complete from previous deploys.

### 5. Create the URL Rewrite Transform Rules
Go to **Rules → Overview → Create rule → URL Rewrite Rule**.

**Rule 1: /periodic**
- Rule name: `Rewrite /periodic to /FlexAIDdS/periodic`
- If incoming requests match: **Custom filter expression**
  ```
  http.host eq "thebonhomme.com" and starts_with(http.request.uri.path, "/periodic")
  ```
- Then:
  - Rewrite to → **Dynamic**
  - Path expression: `concat("/FlexAIDdS", http.request.uri.path)`
  - Query: Preserve query string (checked)
- Deploy

**Rule 2: /entropy-driven** (identical pattern)
- Filter: `... and starts_with(http.request.uri.path, "/entropy-driven")`
- Path: `concat("/FlexAIDdS", http.request.uri.path)`

**Rule 3 (optional but recommended): /entropy-help**
- Same pattern for `/entropy-help`

You can also create a broader rule or use the Wildcard pattern UI if you prefer the visual editor (Request URL `https://thebonhomme.com/periodic/*` → rewrite target `/FlexAIDdS/periodic/${1}`).

### 6. Create the Redirect Rule for /flexaid
This fulfills the original requirement that `/flexaid` "be a redirect (not a page)".

Go to **Rules → Create rule → Redirect Rule** (or "URL Redirect").

- Rule name: `Redirect /flexaid to pre-redesign archive`
- If incoming requests match: Custom filter expression
  ```
  http.host eq "thebonhomme.com" and (http.request.uri.path eq "/flexaid" or starts_with(http.request.uri.path, "/flexaid/"))
  ```
- Then:
  - URL redirect
  - Target URL: `https://thebonhomme.com/archive/pre-redesign-light-dark-2026-05-23/`
  - Status code: `301` (Permanent redirect) — or `302` if you prefer temporary
  - Preserve query string: (optional, not needed here)
- Deploy

The existing `site/flexaid/index.html` (meta refresh + JS redirect) will remain as a fallback but the Cloudflare rule will take precedence and is cleaner.

### 7. Wait for propagation and verify
- Use https://dnschecker.org or run `dig thebonhomme.com` + `dig www.thebonhomme.com` from your machine.
- Full propagation can take up to an hour in rare cases, but is often much faster with Cloudflare.

**Verification commands** (run these):
```bash
curl -I https://thebonhomme.com/periodic
curl -I https://thebonhomme.com/periodic?element=Pt
curl -I https://thebonhomme.com/flexaid
curl -I https://thebonhomme.com/entropy-driven
curl -I https://thebonhomme.com/entropy-help
```

Expected results:
- `/periodic*` → HTTP 200, serves the full interactive table (all the FDA + clinical trial drug data you asked for).
- `/flexaid` → 301/302 redirecting to the archive URL.
- Browser (incognito recommended): address bar must stay at the clean path; content renders correctly (design system, JS search/filters/modals all working).

The old paths (`/FlexAIDdS/periodic` etc.) will continue to work as a direct fallback.

## Post-setup tips
- No action needed in the GitHub Actions workflow for basic functionality.
- If you later want instant cache invalidation on deploys, you can add a Cloudflare API purge step (using a token with Zone.Cache Purge permission) — we can implement that in a follow-up PR once everything is live.
- The site/periodic/README.md and site/flexaid/README.md already document the intended clean URLs.
- Future work (new elements in the table, content updates, etc.) is done exactly as before: edit under `site/`, open PR from a feature branch off `origin/master`, merge → workflow runs → live at clean paths.

## References
- This file was created as the persistent, committed reminder after the repo-side preparation (gitignore + note cleanup) was completed and pushed.
- Original diagnosis came from the raw GitHub Actions log you pasted (the `/FlexAIDdS/` prefix).
- Cloudflare + GitHub Pages is a very common and well-supported pattern.

Once you complete steps 1–7 above, reply here with the output of the `curl -I` commands (or screenshots of the browser) and we will verify together and close the loop on making everything functional at the exact URLs you asked for.

Le Bonhomme Pharma · 2026
