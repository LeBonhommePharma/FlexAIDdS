# B1.1: Register entropy.help Domain + Configure DNS

**Status**: Historical setup plan; current domain availability, ownership, DNS, and deployment are unverified in this repository.

**Recommended Path**:
1. Register on Namecheap (fastest + possible first-year discount)
2. Immediately transfer to Cloudflare Registrar (best long-term pricing, no markup)
3. Point DNS to GitHub Pages (or Vercel) for the static site in `site/entropy-help/`

---

## Step 1: Register the Domain

### Option A (Recommended for speed): Namecheap

1. Go to: https://www.namecheap.com/domains/
2. Search for `entropy.help`
3. Add to cart (first year often heavily discounted for new gTLDs)
4. Complete checkout with privacy protection enabled (free on Namecheap)
5. Note the exact expiry date and renewal price (typically $30–45/year for .help)

### Option B: Cloudflare Registrar (if supported)

1. Log in to Cloudflare dashboard
2. Go to Registrar → Register a domain
3. Search `entropy.help`
4. If available, register at true registry cost (no markup, excellent for long-term)

**Note**: If Cloudflare does not list .help at the time you check, register on Namecheap first, then transfer after 60 days (standard ICANN rule).

---

## Step 2: Set Up DNS (Cloudflare Recommended)

After registration:

### If using Cloudflare Registrar + Cloudflare DNS (best practice)

1. In Cloudflare, add the domain as a site (even if not yet pointed).
2. Cloudflare will give you two nameservers. Change the nameservers at your registrar to Cloudflare's.
3. Once active, create the following records:

**Apex domain (`entropy.help`)**:
- Type: `A`
- Name: `@`
- Value: `185.199.108.153` (GitHub Pages)
- Proxy status: **DNS only** (grey cloud) — required for GitHub Pages custom domains

(You may also add the other GitHub Pages IPs for redundancy:
185.199.109.153
185.199.110.153
185.199.111.153)

**WWW subdomain**:
- Type: `CNAME`
- Name: `www`
- Value: `yourusername.github.io`  (or the actual repo if using project pages)
- Proxy: DNS only

### GitHub Pages Custom Domain Verification

1. Go to your repo → Settings → Pages
2. Under "Custom domain", enter `entropy.help`
3. GitHub will create a `CNAME` file automatically in the `gh-pages` branch (or you can create it manually in the `site/` folder if using a different hosting strategy).

**Important**: Create a file at the root of your published site:
`site/CNAME` containing exactly:
```
entropy.help
```

This tells GitHub Pages to serve the custom domain.

---

## Step 3: HTTPS (Automatic)

- GitHub Pages + custom domain → GitHub automatically provisions a Let's Encrypt certificate (usually within minutes to a few hours).
- Cloudflare can also terminate SSL if you want (but GitHub's is sufficient and simpler).

---

## Step 4: Update References in the Project

Once the domain is live:

1. Update `site/entropy-help/index.html`, `ledger.html`, and `request.html`:
   - Change any `github.io` or relative links to `https://entropy.help/...` where appropriate.
   - Update meta tags and canonical URLs.

2. Update `docs/entropy-help/MANIFESTO.md` and other docs with the final domain.

3. Update the coordination issue and future announcements.

---

## Step 5: Recommended Final Architecture

- **Hosting**: GitHub Pages (free, reliable, automatic deploys from the `site/` folder via GitHub Actions or Pages settings)
- **DNS + CDN**: Cloudflare (best security, DNSSEC, DDoS protection, and cheap/zero markup renewal)
- **Domain owner**: Personal or NRGlab-related entity (document clearly for the project)

---

## Quick Checklist

- [ ] Search and register `entropy.help` on Namecheap
- [ ] Enable WHOIS privacy
- [ ] Add domain to Cloudflare
- [ ] Update nameservers at registrar
- [ ] Create required A + CNAME records (DNS only)
- [ ] Add `CNAME` file in the published site root
- [ ] Verify in GitHub Pages settings
- [ ] Wait for HTTPS certificate
- [ ] Test https://entropy.help and https://www.entropy.help
- [ ] Update all project references
- [ ] Announce the new canonical domain

---

**Historical note**: A May 2026 planning pass described the domain as available, but no registry receipt is deposited here. Re-check availability before acting.

Only a registrar receipt plus live DNS/HTTPS checks can establish completion; this document is not that evidence.
