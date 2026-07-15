# Audit: 5a24ebbd2 — Fix: FlexAID --legacy for A/B/B0 pilot; local-first OUT + deferred iCloud sync

> Synthesized from swarm agent completion for orchestrator (full agent report may live on docs/audit-* branches).

## Summary
FlexAID --legacy + local-first OUT + deferred iCloud sync.

## Severity: CRITICAL

## Findings
### F1. CRITICAL empty ${extra[@]} under set -u on macOS bash 3.2 — launch dies
- Evidence: agent deep audit of `5a24ebbd2`
- Why: science/repro impact
- Fix: see swarm SUMMARY
### F2. Sync arg parser broken (shift corrupts --arm)
- Evidence: agent deep audit of `5a24ebbd2`
- Why: science/repro impact
- Fix: see swarm SUMMARY
### F3. C0 dual-launch guards removed
- Evidence: agent deep audit of `5a24ebbd2`
- Why: science/repro impact
- Fix: see swarm SUMMARY
### F4. RUN_RECEIPT still hard-codes storage icloud_drive while OUT local
- Evidence: agent deep audit of `5a24ebbd2`
- Why: science/repro impact
- Fix: see swarm SUMMARY

## Ranking/scoring impact: NO
## Reproducibility impact: YES
## Tests adequate: NO
## Verdict: BLOCK
