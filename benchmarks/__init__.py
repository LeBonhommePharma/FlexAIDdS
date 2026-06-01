"""Backward-compatibility shim. Real code lives in flexaidds.dataset_runner.

This import is intentionally soft so that subpackages like benchmarks/re-dock/
can be used for thermodynamic / Codex-local benchmark work without requiring
the full flexaidds package to be installed (important on M3 Pro with iCloud-only
storage discipline and partial environments in AI coding sandboxes).
"""

try:
    from flexaidds.dataset_runner import *  # noqa: F401,F403
    from flexaidds.dataset_runner import __all__  # noqa: F401
except ImportError:
    # Allow standalone use of re-dock, nextgen, etc. when flexaidds is not on PYTHONPATH.
    # This is the smallest non-behavioral packaging/import fix needed for the
    # Codex → local iCloud benchmark resume workflow.
    pass
