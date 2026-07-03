from __future__ import annotations

import sys


def _disable_numba_disk_cache() -> None:
    import numba

    def no_disk_cache(decorator):
        def wrapped(*args, **kwargs):
            kwargs["cache"] = False
            return decorator(*args, **kwargs)

        return wrapped

    numba.jit = no_disk_cache(numba.jit)
    numba.njit = no_disk_cache(numba.njit)


def main() -> None:
    _disable_numba_disk_cache()
    from boltz.main import cli

    cli(args=sys.argv[1:])


if __name__ == "__main__":
    main()
