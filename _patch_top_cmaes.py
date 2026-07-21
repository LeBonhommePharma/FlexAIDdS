#!/usr/bin/env python3
# _patch_top_cmaes.py — insert FLEXAIDDS_SEARCH CMA-ES branch into LIB/top.cpp
#
# Invoked by apply_integration.sh. Idempotent: exits 0 if markers already present.
# Branch matches chunk1 adapter API (cmaes_run_dock / fill / write_trace_csv).
#
# Copyright 2026 Le Bonhomme Pharma / Louis-Philippe Morency / NRGlab
# SPDX-License-Identifier: Apache-2.0

from __future__ import annotations

import re
import sys
from pathlib import Path

MARKER_BEGIN = "// FLEXAIDDS_CMAES_BEGIN"
MARKER_END = "// FLEXAIDDS_CMAES_END"

# Indentation matches the live top.cpp "Standard single GA run" block (3 tabs).
BRANCH = r'''			// FLEXAIDDS_CMAES_BEGIN
			// Opt-in CMA-ES search backend. Env FLEXAIDDS_SEARCH=cmaes (or CMAES)
			// swaps the operator; scoring path (ic2cf / 5 seam fns) is unchanged.
			// Eval budget: λ×gens ≡ pop×gens (claim 1000×2000 = 2e6). See CMAES_INTEGRATION.md.
			// Allocation contract matches GA(): sets chrom, chrom_snapshot, gene_lim,
			// memchrom so existing top.cpp free / ranking loops remain valid.
			// ic2cf.cpp / gaboom.cpp are never modified by this branch.
			const char* flexaidds_search = std::getenv("FLEXAIDDS_SEARCH");
			const bool use_cmaes = flexaidds_search &&
				(std::strcmp(flexaidds_search, "cmaes") == 0 ||
				 std::strcmp(flexaidds_search, "CMAES") == 0);
			if (use_cmaes) {
				// Minimal GA-input plumbing so GB->num_chrom / max_generations / seed
				// match the claim budget path (read_gainputs when a .inp is present).
				GB->num_genes = FA->npar;
				if (GB->num_genes <= 0) {
					fprintf(stderr, "ERROR: CMA-ES: no parameters to optimize (FA->npar=0).\n");
					Terminate(1);
				}
				if (gainp[0] != 0) {
					int geninterval = 0, popszpartition = 0;
					read_gainputs(FA, GB, &geninterval, &popszpartition, gainp);
				}
				if (GB->num_chrom <= 0) GB->num_chrom = 1000;
				if (GB->max_generations <= 0) GB->max_generations = 2000;

				gene_lim = (genlim*)malloc(static_cast<size_t>(GB->num_genes) * sizeof(genlim));
				if (!gene_lim) {
					fprintf(stderr, "ERROR: CMA-ES: gene_lim allocation failed.\n");
					Terminate(1);
				}

				CmaesConfig cma_cfg;
				cma_cfg.population = GB->num_chrom > 0 ? GB->num_chrom : 1000;
				cma_cfg.max_evals = 0;
				if (const char* me = std::getenv("FLEXAIDDS_CMAES_MAX_EVALS")) {
					if (me[0] != 0) {
						char* endp = nullptr;
						const long long v = std::strtoll(me, &endp, 10);
						if (endp != me && v > 0)
							cma_cfg.max_evals = static_cast<std::int64_t>(v);
					}
				}
				if (cma_cfg.max_evals <= 0) {
					const long long budget =
						static_cast<long long>(cma_cfg.population) *
						static_cast<long long>(GB->max_generations);
					cma_cfg.max_evals = budget > 0 ? budget : 2000000LL;
				}
				cma_cfg.seed = static_cast<std::uint32_t>(GB->seed != 0 ? GB->seed : 1);
				cma_cfg.enable_entropy_trace = true;
				if (!output_prefix.empty()) {
					cma_cfg.write_trace = output_prefix + "_cmaes_entropy.csv";
				}

				CmaesResult cma_res;
				std::vector<EntropyTraceSample> cma_trace;
				const int cma_rc = cmaes_run_dock(
					FA, GB, VC, gene_lim, atoms, residue, cleftgrid,
					ic2cf, cma_cfg, &cma_res, &cma_trace);
				if (cma_rc != 0) {
					fprintf(stderr, "ERROR: cmaes_run_dock failed (rc=%d status=%d)\n",
					        cma_rc, cma_res.status);
					n_chrom_snapshot = 0;
				} else {
					if (!cma_cfg.write_trace.empty() && !cma_trace.empty()) {
						cmaes_write_trace_csv(cma_cfg.write_trace, cma_trace);
					}
					// Archive size K — each chrom[i].genes / chrom_snapshot[i].genes
					// is a separate malloc so top.cpp free loops stay correct.
					const int arch_n = static_cast<int>(cma_res.archive_genes.size());
					const int K = std::max(1, std::min(cma_cfg.archive_size,
						arch_n > 0 ? arch_n : 1));
					memchrom = K;
					chrom = (chromosome*)calloc(static_cast<size_t>(K), sizeof(chromosome));
					chrom_snapshot = (chromosome*)calloc(static_cast<size_t>(K), sizeof(chromosome));
					// Temporary contiguous storage for cmaes_fill_chromosomes only;
					// immediately re-homed into per-slot malloc buffers.
					gene* tmp_storage = (gene*)calloc(
						static_cast<size_t>(K) * static_cast<size_t>(GB->num_genes),
						sizeof(gene));
					if (!chrom || !chrom_snapshot || !tmp_storage) {
						fprintf(stderr, "ERROR: CMA-ES: chromosome/snapshot allocation failed.\n");
						n_chrom_snapshot = 0;
					} else {
						const int filled = cmaes_fill_chromosomes(
							cma_res, GB->num_genes, chrom, K, tmp_storage);
						n_chrom_snapshot = 0;
						for (int i = 0; i < filled; ++i) {
							gene* cgenes = (gene*)malloc(
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							gene* sgenes = (gene*)malloc(
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							if (!cgenes || !sgenes) {
								fprintf(stderr, "ERROR: CMA-ES: gene buffer allocation failed.\n");
								free(cgenes);
								free(sgenes);
								break;
							}
							std::memcpy(cgenes, chrom[i].genes,
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							std::memcpy(sgenes, chrom[i].genes,
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							chrom[i].genes = cgenes;
							chrom_snapshot[i] = chrom[i];
							chrom_snapshot[i].genes = sgenes;
							++n_chrom_snapshot;
						}
						free(tmp_storage);
						// Patch counters so free(chrom_snapshot) walks K slots, not
						// num_chrom * max_generations (ParallelDock uses the same fix).
						GB->num_chrom = std::max(1, n_chrom_snapshot);
						GB->max_generations = 1;
						printf("[SEARCH] backend=cmaes evals=%d n_snap=%d best_cf=%.6f "
						       "lambda=%d max_evals=%lld\n",
						       cma_res.n_evals, n_chrom_snapshot, cma_res.best_cf,
						       cma_cfg.population,
						       static_cast<long long>(cma_cfg.max_evals));
						if (!cma_cfg.write_trace.empty()) {
							printf("[SEARCH] cmaes entropy trace: %s\n",
							       cma_cfg.write_trace.c_str());
						}
					}
				}
			} else {
				GAContext ga_ctx;
				n_chrom_snapshot = GA(FA,GB,VC,&chrom,&chrom_snapshot,&gene_lim,atoms,residue,&cleftgrid,gainp,&memchrom,ic2cf, &ga_ctx);
			}
			// FLEXAIDDS_CMAES_END
'''


def patch(top_path: Path) -> None:
    text = top_path.read_text(encoding="utf-8")
    if MARKER_BEGIN in text and MARKER_END in text:
        return

    # Primary: whole "Standard single GA run" block.
    # IMPORTANT: use a callable replacement — BRANCH contains '\0' C escapes;
    # a plain re.sub replacement string would expand \0/\1 as group refs.
    pattern = re.compile(
        r"(\t\t\} else \{\n)"
        r"(\t\t\t// ── Standard single GA run ──\n)"
        r"(\t\t\tGAContext ga_ctx;\n)"
        r"(\t\t\tn_chrom_snapshot = GA\(FA,GB,VC,&chrom,&chrom_snapshot,&gene_lim,"
        r"atoms,residue,&cleftgrid,gainp,&memchrom,ic2cf, &ga_ctx\);\n)"
        r"(\t\t\})",
        re.MULTILINE,
    )

    def _repl_primary(m: re.Match[str]) -> str:
        return (
            m.group(1)
            + "\t\t\t// ── Standard single search run (GA default; CMA-ES opt-in) ──\n"
            + BRANCH
            + m.group(5)
        )

    new_text, n = pattern.subn(_repl_primary, text, count=1)

    if n != 1:
        # Fallback: last bare GA call in the standard path.
        pattern2 = re.compile(
            r"^(\t\t\t)n_chrom_snapshot = GA\(FA,GB,VC,&chrom,&chrom_snapshot,"
            r"&gene_lim,atoms,residue,&cleftgrid,gainp,&memchrom,ic2cf, &ga_ctx\);\s*$",
            re.MULTILINE,
        )
        matches = list(pattern2.finditer(text))
        if not matches:
            sys.stderr.write(
                "apply_integration: could not locate standard GA() call in top.cpp\n"
            )
            sys.exit(2)
        m = matches[-1]
        start = m.start()
        # Drop preceding GAContext + comment if present.
        for rx in (
            re.compile(r"\t\t\tGAContext ga_ctx;\s*\n$"),
            re.compile(r"\t\t\t// ── Standard single GA run ──\s*\n$"),
        ):
            prefix = text[max(0, start - 120) : start]
            mm = rx.search(prefix)
            if mm:
                start = start - len(mm.group(0))
        new_text = (
            text[:start]
            + "\t\t\t// ── Standard single search run (GA default; CMA-ES opt-in) ──\n"
            + BRANCH
            + text[m.end() :]
        )

    if MARKER_BEGIN not in new_text or MARKER_END not in new_text:
        sys.stderr.write("apply_integration: markers missing after patch\n")
        sys.exit(3)

    top_path.write_text(new_text, encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(f"usage: {argv[0]} LIB/top.cpp", file=sys.stderr)
        return 2
    path = Path(argv[1])
    if not path.is_file():
        print(f"missing {path}", file=sys.stderr)
        return 2
    patch(path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
