# Run monitoring — astex_repro/full
- Worker: benchmark_datasets (reparented, survives bash-wrapper death)
- LIVENESS: stderr mtime advancing, NOT the wrapper PID in run.pid
  `find full -name stderr.log -newermt '-3 min' | head` → nonempty = alive
- PROGRESS: `ls full/*/result.csv | wc -l` → N/85 targets done
- RESUME after any death: `./run_full.sh` (no --force → skips done targets)
- Engine PINNED: engine/FlexAIDdS (sha 9f47c1bb) + staged .dat/.def data files
- SCORING: python3 score_reference.py → poster_metric_reference.csv (spyrmsd, PoseBusters-grade symmetry RMSD, pooled <2A)
