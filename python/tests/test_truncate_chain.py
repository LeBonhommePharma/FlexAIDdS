"""Tests for flexaidds.truncate_chain (nascent-chain PDB writer)."""

from __future__ import annotations

from pathlib import Path

import pytest

from flexaidds.truncate_chain import write_extended_ca_chain, main, _CA_CA_DIST_A


def test_write_extended_ca_chain_geometry(tmp_path: Path):
    out = tmp_path / "chain.pdb"
    seq = "ACDEFGHIKLMNPQRSTVWY"
    result = write_extended_ca_chain(seq, length=5, output_path=out)
    assert result == out
    text = out.read_text()
    assert "Synthetic nascent chain at L_k = 5" in text
    atom_lines = [ln for ln in text.splitlines() if ln.startswith("ATOM")]
    assert len(atom_lines) == 5
    # Residues: A C D E F → ALA CYS ASP GLU PHE
    assert "ALA" in atom_lines[0]
    assert "PHE" in atom_lines[4]
    # Cα spacing along +X
    xs = []
    for ln in atom_lines:
        xs.append(float(ln[30:38]))
    for i in range(1, len(xs)):
        assert xs[i] - xs[i - 1] == pytest.approx(_CA_CA_DIST_A)
    assert text.strip().endswith("END")


def test_write_truncates_to_length(tmp_path: Path):
    out = tmp_path / "short.pdb"
    write_extended_ca_chain("AAAAAA", length=2, output_path=out)
    atoms = [ln for ln in out.read_text().splitlines() if ln.startswith("ATOM")]
    assert len(atoms) == 2


def test_unknown_residue_becomes_unk(tmp_path: Path):
    out = tmp_path / "unk.pdb"
    write_extended_ca_chain("AZ", length=2, output_path=out)
    text = out.read_text()
    assert "UNK" in text  # Z is not a standard AA code in the map


def test_cli_writes_pdb(tmp_path: Path):
    out = tmp_path / "cli.pdb"
    rc = main(["--sequence", "GAGA", "--length", "3", "--output", str(out)])
    assert rc == 0
    assert out.is_file()
    atoms = [ln for ln in out.read_text().splitlines() if ln.startswith("ATOM")]
    assert len(atoms) == 3


def test_cli_reads_fasta(tmp_path: Path):
    fa = tmp_path / "seq.fa"
    fa.write_text(">chain\nMKTAY\nIAK\n")
    out = tmp_path / "from_fa.pdb"
    rc = main(["--sequence", f"@{fa}", "--length", "6", "--output", str(out)])
    assert rc == 0
    atoms = [ln for ln in out.read_text().splitlines() if ln.startswith("ATOM")]
    assert len(atoms) == 6


def test_cli_rejects_empty_and_bad_length(tmp_path: Path):
    out = tmp_path / "bad.pdb"
    assert main(["--sequence", "", "--length", "1", "--output", str(out)]) == 2
    assert main(["--sequence", "AAA", "--length", "0", "--output", str(out)]) == 2
    assert main(["--sequence", "AAA", "--length", "9", "--output", str(out)]) == 2
    assert not out.exists()
