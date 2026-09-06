"""Tests for flexaidds.__main__ – CLI entry point.

Covers build_parser() and the output branches of main():
  - human-readable summary (default)
  - machine-readable JSON (--json flag)
  - CSV export (--csv flag)
  - top-N filtering (--top flag)
  - version display (--version flag)
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

from flexaidds.__main__ import build_parser, main
from flexaidds.__version__ import __version__
from flexaidds.models import DockingResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_pdb(path: Path, remarks: list[str]) -> None:
    lines = [f"REMARK {r}\n" for r in remarks]
    lines += [
        "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n",
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


# ===========================================================================
# build_parser
# ===========================================================================

class TestBuildParser:
    def test_returns_argument_parser(self):
        import argparse
        parser = build_parser()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_results_dir_positional(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([str(tmp_path)])
        assert args.results_dir == tmp_path

    def test_json_flag_default_false(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([str(tmp_path)])
        assert args.json is False

    def test_json_flag_enabled(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([str(tmp_path), "--json"])
        assert args.json is True

    def test_prog_name(self):
        parser = build_parser()
        assert "flexaidds" in parser.prog

    def test_version_flag(self, capsys):
        parser = build_parser()
        with pytest.raises(SystemExit, match="0"):
            parser.parse_args(["--version"])
        out = capsys.readouterr().out
        from flexaidds.__version__ import __version__
        assert __version__ in out


# ===========================================================================
# main() – human-readable output
# ===========================================================================

class TestMainHumanOutput:
    def _make_dir(self, tmp_path: Path) -> Path:
        _write_pdb(
            tmp_path / "mode_1_pose_1.pdb",
            ["binding_mode = 1", "pose_rank = 1", "CF = -42.5",
             "free_energy = -41.0", "temperature = 300.0"],
        )
        _write_pdb(
            tmp_path / "mode_1_pose_2.pdb",
            ["binding_mode = 1", "pose_rank = 2", "CF = -39.0",
             "temperature = 300.0"],
        )
        _write_pdb(
            tmp_path / "mode_2_pose_1.pdb",
            ["binding_mode = 2", "pose_rank = 1", "CF = -35.0",
             "temperature = 300.0"],
        )
        return tmp_path

    def test_returns_zero(self, tmp_path, monkeypatch):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d)])
        assert main() == 0

    def test_prints_n_modes(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d)])
        main()
        out = capsys.readouterr().out
        assert "2" in out  # n_modes

    def test_prints_temperature(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d)])
        main()
        out = capsys.readouterr().out
        assert "300" in out

    def test_prints_results_directory(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d)])
        main()
        out = capsys.readouterr().out
        assert str(d) in out

    def test_prints_top_mode_info(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d)])
        main()
        out = capsys.readouterr().out
        assert "mode_id" in out

    def test_no_temperature_line_when_absent(self, tmp_path, monkeypatch, capsys):
        _write_pdb(tmp_path / "mode_1_pose_1.pdb", ["CF = -5.0"])
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(tmp_path)])
        main()
        out = capsys.readouterr().out
        assert "Temperature" not in out

    def test_no_top_mode_line_when_no_modes(self, tmp_path, monkeypatch, capsys):
        """Covered via empty-directory – load_results raises; tested separately."""
        # With a single unscored pose top_mode() still finds a mode, so just
        # verify the "Top mode" label appears only when top_mode() returns non-None.
        _write_pdb(tmp_path / "mode_1_pose_1.pdb", ["CF = -5.0"])
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(tmp_path)])
        main()
        out = capsys.readouterr().out
        assert "Top mode" in out


# ===========================================================================
# main() – JSON output
# ===========================================================================

class TestMainJsonOutput:
    def _make_dir(self, tmp_path: Path) -> Path:
        _write_pdb(
            tmp_path / "mode_1_pose_1.pdb",
            ["binding_mode = 1", "pose_rank = 1", "CF = -42.5",
             "free_energy = -41.0", "temperature = 300.0"],
        )
        return tmp_path

    def test_returns_zero(self, tmp_path, monkeypatch):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--json"])
        assert main() == 0

    def test_output_is_valid_json(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--json"])
        main()
        out = capsys.readouterr().out
        parsed = json.loads(out)  # raises if not valid JSON
        assert isinstance(parsed, dict)

    def test_json_has_required_keys(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--json"])
        main()
        out = capsys.readouterr().out
        parsed = json.loads(out)
        assert set(parsed.keys()) >= {"source_dir", "n_modes", "binding_modes",
                                       "temperature", "metadata"}

    def test_json_n_modes(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--json"])
        main()
        parsed = json.loads(capsys.readouterr().out)
        assert parsed["n_modes"] == 1

    def test_json_temperature(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--json"])
        main()
        parsed = json.loads(capsys.readouterr().out)
        assert parsed["temperature"] == pytest.approx(300.0)

    def test_json_binding_modes_list(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--json"])
        main()
        parsed = json.loads(capsys.readouterr().out)
        assert isinstance(parsed["binding_modes"], list)
        assert len(parsed["binding_modes"]) == 1

    def test_json_source_dir_is_string(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--json"])
        main()
        parsed = json.loads(capsys.readouterr().out)
        assert isinstance(parsed["source_dir"], str)


# ===========================================================================
# --version flag
# ===========================================================================

class TestVersionFlag:
    def test_version_flag_exits_zero(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["flexaidds", "--version"])
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 0

    def test_version_flag_prints_version(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["flexaidds", "--version"])
        with pytest.raises(SystemExit):
            main()
        out = capsys.readouterr().out
        assert __version__ in out

    def test_parser_accepts_version(self):
        parser = build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["--version"])
        assert exc_info.value.code == 0

    def test_short_version_flag(self, capsys):
        parser = build_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args(["-V"])
        assert exc_info.value.code == 0


# ===========================================================================
# main() – CSV output
# ===========================================================================

class TestMainCsvOutput:
    def _make_dir(self, tmp_path: Path) -> Path:
        _write_pdb(
            tmp_path / "mode_1_pose_1.pdb",
            ["binding_mode = 1", "pose_rank = 1", "CF = -42.5",
             "free_energy = -41.0", "temperature = 300.0"],
        )
        _write_pdb(
            tmp_path / "mode_2_pose_1.pdb",
            ["binding_mode = 2", "pose_rank = 1", "CF = -35.0",
             "temperature = 300.0"],
        )
        return tmp_path

    def test_returns_zero(self, tmp_path, monkeypatch):
        d = self._make_dir(tmp_path)
        csv_path = tmp_path / "out.csv"
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--csv", str(csv_path)])
        assert main() == 0

    def test_writes_csv_file(self, tmp_path, monkeypatch):
        d = self._make_dir(tmp_path)
        csv_path = tmp_path / "out.csv"
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--csv", str(csv_path)])
        main()
        assert csv_path.exists()

    def test_csv_has_header_and_rows(self, tmp_path, monkeypatch):
        d = self._make_dir(tmp_path)
        csv_path = tmp_path / "out.csv"
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--csv", str(csv_path)])
        main()
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
        assert len(rows) == 2
        assert "mode_id" in rows[0]
        assert "free_energy" in rows[0]

    def test_csv_round_trips_via_from_csv(self, tmp_path, monkeypatch):
        d = self._make_dir(tmp_path)
        csv_path = tmp_path / "out.csv"
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--csv", str(csv_path)])
        main()

        restored = DockingResult.from_csv(csv_path)
        assert restored.n_modes == 2
        assert restored.binding_modes[0].mode_id == 1
        assert restored.binding_modes[0].free_energy == pytest.approx(-41.0)

    def test_prints_confirmation(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        csv_path = tmp_path / "out.csv"
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--csv", str(csv_path)])
        main()
        out = capsys.readouterr().out
        assert "2" in out
        assert str(csv_path) in out


# ===========================================================================
# main() – --top flag
# ===========================================================================

class TestMainTopFlag:
    def _make_dir(self, tmp_path: Path) -> Path:
        for i in range(1, 4):
            _write_pdb(
                tmp_path / f"mode_{i}_pose_1.pdb",
                [f"binding_mode = {i}", "pose_rank = 1",
                 f"CF = {-40.0 + i}", f"free_energy = {-39.0 + i}",
                 "temperature = 300.0"],
            )
        return tmp_path

    def test_top_limits_table_rows(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d), "--top", "1"])
        main()
        out = capsys.readouterr().out

        # Header reports the TRUE total, not the truncated one. Colour
        # accessors return "" off a TTY, so under capsys this is
        # "binding modes=3" (26f8ee17 design-system kv form).
        assert "binding modes=3" in out

        # --top 1 must leave exactly one data row.
        data_rows = [ln for ln in out.splitlines() if ln.strip()[:1].isdigit()]
        assert len(data_rows) == 1, (
            f"--top 1 should leave exactly 1 data row, got {len(data_rows)}: {data_rows}"
        )

    def test_top_none_shows_all_rows(self, tmp_path, monkeypatch, capsys):
        d = self._make_dir(tmp_path)
        monkeypatch.setattr(sys, "argv", ["flexaidds", str(d)])
        main()
        out = capsys.readouterr().out
        assert "binding modes=3" in out
        data_rows = [ln for ln in out.splitlines() if ln.strip()[:1].isdigit()]
        assert len(data_rows) == 3, (
            f"no --top should list all 3 modes, got {len(data_rows)}: {data_rows}"
        )

    def test_top_default_is_none(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([str(tmp_path)])
        assert args.top is None

    def test_top_parsed_as_int(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([str(tmp_path), "--top", "5"])
        assert args.top == 5
