# thebonhomme.com — Cloudflare & DNS setup (beginner guide)

**Last updated:** 2026-06-22  
**Read this if:** FlexAID∆S, Entropy Docking, or Mol* look broken, or you see a stub homepage.

---

## What Cloudflare is (30 seconds)

When someone types `thebonhomme.com`, their browser asks DNS “where is this site?” Cloudflare (if you use it) sits between the visitor and GitHub Pages. It can:

- Point the domain to the right server (DNS)
- Speed up / cache pages (CDN — orange cloud)
- Rewrite URLs (Transform Rules) — **this is what broke your site before**

**Today your site files live on GitHub.** One repo owns the domain:  
`LeBonhommePharma/lebonhommepharma.github.io` → **thebonhomme.com**

A second repo (`FlexAIDdS` `gh-pages`) must **not** also claim `thebonhomme.com` or `/FlexAIDdS/` shows the wrong page.

---

## Part A — Fix GitHub first (do this before Cloudflare)

### A1. Confirm who serves thebonhomme.com

1. Open https://github.com/LeBonhommePharma/lebonhommepharma.github.io/settings/pages  
2. **Custom domain** must be: `thebonhomme.com`  
3. **Enforce HTTPS** should be on.

### A2. Remove duplicate domain from FlexAIDdS gh-pages

1. Open https://github.com/LeBonhommePharma/FlexAIDdS/tree/gh-pages  
2. If you see a file named **`CNAME`** containing `thebonhomme.com` → **delete it** (commit to `gh-pages`).  
3. Future deploys: our workflow now excludes `CNAME` from `gh-pages` automatically.

**Why:** Two repos with the same `CNAME` made `/FlexAIDdS/` serve the corporate homepage instead of the product page.

### A3. Full site sync (already automated)

CI runs `scripts/sync_apex_to_usersite.sh` after each deploy. It copies the entire `site/` folder to `lebonhommepharma.github.io`. That repo is what visitors get at `/`, `/FlexAIDdS/`, `/entropy-driven/`, etc.

---

## Part B — Cloudflare (step by step)

### Do you even use Cloudflare?

On your Mac, open **Terminal** and run:

```bash
dig thebonhomme.com +short
```

| You see | Meaning |
|--------|---------|
| `185.199.108.153` (and `.109`, `.110`, `.111`) | DNS goes **straight to GitHub**. Cloudflare proxy rules are **not** active unless you changed nameservers. |
| `104.x` or `172.x` Cloudflare IPs | Domain is **proxied through Cloudflare**. Follow all steps below. |

---

### B1. Log in

1. Go to https://dash.cloudflare.com  
2. Sign in (or create a free account).  
3. Click your site **`thebonhomme.com`**.  
   - If it’s **not listed**, you’re probably on GitHub-only DNS (Part B optional). Skip to **Part C — Verify**.

---

### B2. Delete ALL old Transform Rules (critical)

Old rules prepended `/FlexAIDdS` to every path. That was for an outdated deploy layout and **breaks the site now**.

1. Left sidebar → **Rules** → **Overview** (or **Transform Rules**).  
2. Open **URL Rewrite Rules** (and **Redirect Rules** if listed).  
3. **Delete** every rule whose name or path mentions any of:
   - `Rewrite apex homepage to /FlexAIDdS`
   - `Rewrite /periodic to /FlexAIDdS`
   - `Rewrite /entropy-driven`
   - `concat("/FlexAIDdS"`
   - Redirect `/flexaid` to `archive/pre-redesign`
4. Click **Deploy** / confirm deletes.

**You should have zero URL Rewrite rules** for this site with the current architecture.

---

### B3. DNS records (if Cloudflare manages your domain)

1. Left sidebar → **DNS** → **Records**.  
2. Set this up (orange cloud = Proxied, grey = DNS only):

| Type | Name | Content | Proxy |
|------|------|---------|-------|
| A | `@` | `185.199.108.153` | Your choice* |
| A | `@` | `185.199.109.153` | Your choice* |
| A | `@` | `185.199.110.153` | Your choice* |
| A | `@` | `185.199.111.153` | Your choice* |
| CNAME | `www` | `lebonhommepharma.github.io` | Your choice* |

\* **Grey cloud (DNS only)** = simplest; behaves like today, no rewrite rules needed.  
**Orange cloud (Proxied)** = Cloudflare CDN; still fine **if you deleted all Transform Rules**.

3. Delete duplicate or conflicting A/CNAME records for `@` or `www`.

---

### B4. Optional redirect rules (only these two)

Left sidebar → **Rules** → **Redirect Rules** → **Create rule**.

**Rule 1 — www → apex**

| Field | Value |
|-------|-------|
| Name | `www to apex` |
| Expression | `(http.host eq "www.thebonhomme.com")` |
| Type | Dynamic |
| URL | `concat("https://thebonhomme.com", http.request.uri.path)` |
| Status | `301` |

**Rule 2 — legacy /flexaid shortcut**

| Field | Value |
|-------|-------|
| Name | `flexaid to FlexAIDdS` |
| Expression | `(http.host eq "thebonhomme.com" and http.request.uri.path eq "/flexaid") or (http.host eq "thebonhomme.com" and starts_with(http.request.uri.path, "/flexaid/"))` |
| Target URL | `https://thebonhomme.com/FlexAIDdS/` |
| Status | `301` |

Do **not** add rules that rewrite `/periodic`, `/entropy-driven`, `/drug-of-the-day`, or `/FlexAIDdS` to a different path. GitHub already has the correct folders.

---

### B5. Purge cache (after rule changes)

1. Left sidebar → **Caching** → **Configuration**.  
2. Click **Purge Everything** → confirm.  
3. Wait 2–5 minutes.

---

## Part C — Verify (copy-paste in Terminal)

```bash
# Corporate homepage (full brand page, NOT the stub)
curl -sL https://thebonhomme.com/ | grep -o 'Bonhomme disagrees'

# FlexAID∆S product (React) — title must mention FlexAID∆S, NOT "Le Bonhomme Pharma" only
curl -sL https://thebonhomme.com/FlexAIDdS/ | grep '<title>'

# Entropy-driven landing
curl -sI https://thebonhomme.com/entropy-driven/ | head -1

# Entropy Docking alias → entropy-driven
curl -sI https://thebonhomme.com/entropy-docking/ | head -1

# Drug of the Day
curl -sI https://thebonhomme.com/drug-of-the-day/ | head -1
```

**Pass criteria:**

| URL | Expected |
|-----|----------|
| `/` | Contains `Bonhomme disagrees`; ~33 KB page with manifesto & products |
| `/FlexAIDdS/` | `<title>` contains `FlexAID∆S`; file size ~2–3 KB (not 33 KB) |
| `/entropy-driven/` | HTTP `200` |
| `/entropy-docking/` | HTTP `200` (redirect page) |
| `/drug-of-the-day/` | HTTP `200` |

Open in browser (incognito):  
https://thebonhomme.com/  
https://thebonhomme.com/FlexAIDdS/  
https://thebonhomme.com/entropy-driven/

---

## Part D — URL map (no Cloudflare magic required)

| Clean URL | What it is |
|-----------|------------|
| `/` | Le Bonhomme Pharma corporate homepage |
| `/FlexAIDdS/` | FlexAID∆S product site (React) |
| `/entropy-driven/` | Entropy-driven docking landing |
| `/entropy-docking/` | Alias → `/entropy-driven/` |
| `/entropy-help/` | Entropy.help audits |
| `/drug-of-the-day/` | Drug of the Day series |
| `/periodic/` | Periodic table of psychoactives |
| `/flexaid/` | Redirect → `/FlexAIDdS/` |

All paths are real folders in `site/` on `lebonhommepharma.github.io`. **No path prefix `/FlexAIDdS/` rewrite is needed.**

---

## Part E — Mol* 3D background

The Mol* viewer lives on the **legacy apex** bundle (`site/app.js` + `#molstar-viewer`). The corporate homepage at `/` intentionally does **not** include Mol* (red brand shell). For the molecular hero:

- Use **`/FlexAIDdS/`** (product) or the archived flexaid landing once Mol* is wired into the React hero.

---

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| Stub “site update in progress” | User-site out of date → re-run deploy or `GITHUB_TOKEN=$(gh auth token) bash scripts/sync_apex_to_usersite.sh` |
| `/FlexAIDdS/` shows corporate homepage | Delete `CNAME` on FlexAIDdS `gh-pages`; purge Cloudflare cache |
| CI deploy fails on “stats commit” | Non-fatal now; deploy still continues. Re-run workflow if needed. |
| Old paths like `/FlexAIDdS/periodic` | Use `/periodic/` instead |

---

## Nameservers (only if you moved DNS to Cloudflare)

If your registrar still points to Cloudflare nameservers, keep them. If you never moved off the registrar’s default nameservers and `dig` shows `185.199.x`, you’re on **GitHub Pages DNS** — Parts B2–B5 apply only if you later add the domain to Cloudflare.

---

Le Bonhomme Pharma · Montréal · 2026