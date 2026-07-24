// Audit every local benchmark/canary version: success rates + code modifications
// between git pins. Local-first only (never CloudDocs tree walks).
//
// Invoke:
//   /audit-benchmark-versions
//   /workflow audit-benchmark-versions
//   /workflow audit-benchmark-versions {"roots":"~/flexaidds_results","out_dir":"~/flexaidds_results/workorders"}
//
// Deterministic inventory is scripts/audit_benchmark_versions.py; agents interpret
// and flag regressions. Evidence gate is fail-closed on missing inventory output.

let meta = #{
    name: "audit-benchmark-versions",
    description: "Inventory all local benchmark/canary versions, success rates (genuine/S1), and code changes between git pins",
    when_to_use: "When you need a cross-version audit of docking benchmarks, canary success rates, matrix pins, and which commits changed between runs",
    phases: [
        #{ title: "Inventory", detail: "Run audit_benchmark_versions.py over local result roots" },
        #{ title: "Interpret", detail: "Summarize rates, pin deltas, and regressions" },
        #{ title: "Verify", detail: "Cross-check table counts against JSON" },
    ],
};

let inventory_schema = #{
    "type": "object",
    "required": ["ok", "json_path", "md_path", "n_campaigns", "summary"],
    "properties": #{
        "ok": #{ "type": "boolean" },
        "json_path": #{ "type": "string" },
        "md_path": #{ "type": "string" },
        "n_campaigns": #{ "type": "integer" },
        "summary": #{ "type": "string" },
        "error": #{ "type": "string" },
    },
};

let interpret_schema = #{
    "type": "object",
    "required": ["headline", "best_rate", "worst_regressions", "key_modifications", "gate_status"],
    "properties": #{
        "headline": #{ "type": "string" },
        "best_rate": #{ "type": "string" },
        "worst_regressions": #{
            "type": "array",
            "maxItems": 12,
            "items": #{ "type": "string" },
        },
        "key_modifications": #{
            "type": "array",
            "maxItems": 20,
            "items": #{ "type": "string" },
        },
        "gate_status": #{ "type": "string" },
        "notes": #{ "type": "string" },
    },
};

let verify_schema = #{
    "type": "object",
    "required": ["real", "reason", "evidence"],
    "properties": #{
        "real": #{ "type": "boolean" },
        "reason": #{ "type": "string" },
        "evidence": #{ "type": "string" },
    },
};

// ── args (optional overrides) ──────────────────────────────────────────────
let roots = if args == () { () } else { args.roots };
let out_dir = if args == () { () } else { args.out_dir };
let repo = if args == () { () } else { args.repo };

if roots == () { roots = "~/flexaidds_results"; }
if out_dir == () { out_dir = "~/flexaidds_results/workorders"; }
if repo == () { repo = "."; }

// ── Phase 1: deterministic inventory ───────────────────────────────────────
phase("Inventory");
log("Running scripts/audit_benchmark_versions.py (local roots only)");

let inv_prompt = "";
inv_prompt += "You are running a deterministic inventory of FlexAIDdS benchmark versions.\n";
inv_prompt += "Repository root for git correlation: " + repo + "\n";
inv_prompt += "Result roots (local APFS only — NEVER walk iCloud/Mobile Documents): " + roots + "\n";
inv_prompt += "Output directory: " + out_dir + "\n\n";
inv_prompt += "REQUIRED steps (use the shell / run_terminal_command tool):\n";
inv_prompt += "1. cd to the FlexAIDdS git repo (find AGENTS.md + scripts/audit_benchmark_versions.py).\n";
inv_prompt += "2. Run exactly:\n";
inv_prompt += "   python3 scripts/audit_benchmark_versions.py --roots " + roots;
inv_prompt += " --out-dir " + out_dir + " --repo " + repo + "\n";
inv_prompt += "3. Confirm these files exist:\n";
inv_prompt += "   - <out-dir>/audit_benchmark_versions.json\n";
inv_prompt += "   - <out-dir>/audit_benchmark_versions.md\n";
inv_prompt += "4. Read the printed summary lines and the JSON n_campaigns field.\n\n";
inv_prompt += "Return JSON: ok=true only if both files were written and the script exit code was 0.\n";
inv_prompt += "json_path and md_path must be absolute paths. summary = first ~30 lines of the MD table or CLI summary.\n";
inv_prompt += "If the script fails, ok=false and put stderr in error. Do not invent campaign counts.";

let inv = agent(inv_prompt, #{
    label: "inventory",
    capability_mode: "execute",
    output_schema: inventory_schema,
});

// Fail closed only on hard failure or explicit ok=false. Missing ok (canned smoke
// host) continues with default paths under out_dir so later phases still compile.
if inv == () || !inv.success || inv.output == () {
    complete(#{
        ok: false,
        stage: "Inventory",
        error: "inventory agent failed or returned empty",
        hint: "Run: python3 scripts/audit_benchmark_versions.py --roots ~/flexaidds_results",
    });
}
if inv.output.ok == false {
    let err = "inventory reported ok=false";
    if inv.output.error != () { err = inv.output.error; }
    complete(#{
        ok: false,
        stage: "Inventory",
        error: err,
        hint: "Run: python3 scripts/audit_benchmark_versions.py --roots ~/flexaidds_results",
    });
}

let json_path = out_dir + "/audit_benchmark_versions.json";
let md_path = out_dir + "/audit_benchmark_versions.md";
let n_camps = 0;
if inv.output.json_path != () && inv.output.json_path != "" {
    json_path = inv.output.json_path;
}
if inv.output.md_path != () && inv.output.md_path != "" {
    md_path = inv.output.md_path;
}
if inv.output.n_campaigns != () {
    n_camps = inv.output.n_campaigns;
}
log("Inventory OK: n_campaigns=" + n_camps.to_string() + " md=" + md_path);

// ── Phase 2: interpret success rates + modifications ───────────────────────
phase("Interpret");

let int_prompt = "";
int_prompt += "Read the benchmark version audit artifacts and produce a concise scientific summary.\n";
int_prompt += "Files (read with read_file — do not invent):\n";
int_prompt += "  MD:   " + md_path + "\n";
int_prompt += "  JSON: " + json_path + "\n\n";
int_prompt += "Rules:\n";
int_prompt += "- genuine = seed_echo=0 + non-seed pose + RMSD≤2Å (claim-adjacent canary metric).\n";
int_prompt += "- S1 = rmsd_to_crystal≤2 only (diagnostic; NEVER call it claim success).\n";
int_prompt += "- Local-first: ignore any iCloud paths.\n";
int_prompt += "- Quote actual rates and commit SHAs from the files.\n";
int_prompt += "- From version_deltas / 'Modifications between benchmark pins', list the most\n";
int_prompt += "  important code changes that coincided with rate moves.\n";
int_prompt += "- worst_regressions: campaigns or pin steps where genuine rate is 0 or dropped.\n";
int_prompt += "- gate_status: whether any canary reached genuine ≥1/3 (PASS/FAIL wording).\n";
int_prompt += "- headline: one sentence overall state of docking power across versions.\n";
int_prompt += "Empty arrays are valid only after reading the files.";

let interp = agent(int_prompt, #{
    label: "interpret",
    capability_mode: "read-only",
    output_schema: interpret_schema,
});

// ── Phase 3: adversarial verify counts ─────────────────────────────────────
phase("Verify");

let ver_prompt = "";
ver_prompt += "Adversarially verify the interpret summary against the inventory JSON.\n";
ver_prompt += "JSON path: " + json_path + "\n";
ver_prompt += "MD path: " + md_path + "\n";
if interp != () && interp.success && interp.output != () {
    ver_prompt += "Claimed headline: " + json_encode(interp.output.headline) + "\n";
    ver_prompt += "Claimed gate_status: " + json_encode(interp.output.gate_status) + "\n";
    ver_prompt += "Claimed best_rate: " + json_encode(interp.output.best_rate) + "\n";
}
ver_prompt += "\nSet real=true ONLY if you independently opened the JSON/MD and confirmed that:\n";
ver_prompt += "(1) n_campaigns matches the table, (2) genuine rates cited are present in the data,\n";
ver_prompt += "(3) gate_status wording is consistent with genuine ≥1/3 on any canary.\n";
ver_prompt += "Otherwise real=false. Put concrete file evidence in evidence.";

let ver = agent(ver_prompt, #{
    label: "verify",
    capability_mode: "read-only",
    output_schema: verify_schema,
});

let verified = false;
let ver_reason = "verify agent missing";
let ver_evidence = "";
if ver != () && ver.success && ver.output != () {
    verified = ver.output.real == true;
    if ver.output.reason != () { ver_reason = ver.output.reason; }
    if ver.output.evidence != () { ver_evidence = ver.output.evidence; }
}

// ── write scratch report + complete ────────────────────────────────────────
let report = "";
report += "# Benchmark version audit — workflow result\n\n";
report += "- Inventory: " + n_camps.to_string() + " campaigns\n";
report += "- JSON: `" + json_path + "`\n";
report += "- MD: `" + md_path + "`\n";
report += "- Verified: " + verified.to_string() + "\n\n";
if interp != () && interp.success && interp.output != () {
    report += "## Headline\n\n" + interp.output.headline + "\n\n";
    report += "## Best rate\n\n" + interp.output.best_rate + "\n\n";
    report += "## Gate\n\n" + interp.output.gate_status + "\n\n";
    report += "## Key modifications\n\n";
    if interp.output.key_modifications != () {
        for m in interp.output.key_modifications {
            report += "- " + m + "\n";
        }
    }
    report += "\n## Worst regressions\n\n";
    if interp.output.worst_regressions != () {
        for r in interp.output.worst_regressions {
            report += "- " + r + "\n";
        }
    }
}
report += "\n## Verification\n\n";
report += "- real: " + verified.to_string() + "\n";
report += "- reason: " + ver_reason + "\n";
report += "- evidence: " + ver_evidence + "\n";

let scratch = write_scratch_file("audit_benchmark_versions_report.md", report);

complete(#{
    ok: true,
    n_campaigns: n_camps,
    json_path: json_path,
    md_path: md_path,
    verified: verified,
    report_scratch: scratch,
    summary: if interp != () && interp.success && interp.output != () {
        interp.output.headline
    } else {
        "Inventory complete; interpret unavailable"
    },
});
