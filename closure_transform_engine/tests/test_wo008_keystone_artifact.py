"""WO-008 keystone artifact tests.

Each test exercises one observable claim from the 24-600 lambda=12
spectral-bridge result. The fixture `result` runs the keystone protocol
once and reuses it; all dimension / residual / representation checks
are then read out of the structured KeystoneResult.

These tests are intentionally tolerant on numerical noise (1e-10) but
strict on structural integers (10, 15, 25) and verdict strings.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]

# The docs/figures/outputs for the spectral bridge can live in three places:
#   - main VFD research repo: docs/keystone_24_600/
#   - combined two-paper public bundle: docs/02-spectral-bridge/
#   - legacy single-paper bundle: docs/
_DOCS_CANDIDATES = [
    REPO / "docs" / "keystone_24_600",
    REPO / "docs" / "02-spectral-bridge",
    REPO / "docs",
]
DOCS = next((d for d in _DOCS_CANDIDATES if (d / "spectral_bridge_note.md").exists()),
            _DOCS_CANDIDATES[0])
OUTPUTS = DOCS / "outputs"

# The README lives at DOCS/README.md (main repo) or at the bundle root.
_README_CANDIDATES = [DOCS / "README.md", REPO / "README.md"]
README_PATH = next((p for p in _README_CANDIDATES if p.exists()), _README_CANDIDATES[0])

if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from closure_transform_engine.keystone import run_keystone


@pytest.fixture(scope="module")
def result():
    return run_keystone()


# --- 1. Runner integration ---------------------------------------------

def test_runner_executes_without_error():
    """The keystone runner script exits 0 and produces the expected artifacts."""
    script = REPO / "closure_transform_engine" / "examples" / "run_wo008_keystone.py"
    proc = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        cwd=str(REPO),
    )
    assert proc.returncode == 0, f"runner failed:\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    assert "SELECTIVE_SPECTRAL_BRIDGE" in proc.stdout


# --- 2. JSON summary ----------------------------------------------------

def test_summary_json_present_and_well_formed():
    p = OUTPUTS / "wo008_summary.json"
    assert p.exists(), f"missing: {p}"
    with p.open() as f:
        data = json.load(f)
    for key in (
        "n_vertices", "coset_sizes", "lambda12_lifted_dim",
        "lambda12_global_dim", "lambda12_residual_dim",
        "lambda12_max_residual", "verdict", "summary_strings",
    ):
        assert key in data, f"summary.json missing key: {key}"
    assert data["verdict"] == "SELECTIVE_SPECTRAL_BRIDGE"


# --- 3. CSV verdict table ----------------------------------------------

def test_verdict_csv_present_and_has_all_rows():
    p = OUTPUTS / "wo008_verdict_table.csv"
    assert p.exists(), f"missing: {p}"
    text = p.read_text()
    # Expect exactly one FULL row (lambda=12), one PARTIAL (lambda=0), and
    # three EMPTY rows (lambda=4,8,10).
    assert "FULL" in text
    assert "PARTIAL" in text
    assert "EMPTY" in text
    rows = [r for r in text.splitlines() if r.strip()]
    assert len(rows) == 6, f"expected exactly 6 lines (header + 5 rows), got {len(rows)}"


# --- 4. lambda=12 residual is below tolerance --------------------------

def test_lambda12_residual_below_threshold(result):
    assert result.lambda12_max_residual < 1e-10, (
        f"lambda=12 lift residual {result.lambda12_max_residual:.3e} not < 1e-10"
    )


# --- 5. negative controls remain above threshold -----------------------

def test_negative_controls_all_above_threshold(result):
    for c in result.negative_controls:
        assert c["max_residual"] > 1.0, (
            f"control lambda={c['lambda']}: residual {c['max_residual']:.3e} "
            f"unexpectedly small (selectivity broken)"
        )


def test_selectivity_via_intersection_dim(result):
    """The lift image intersects E_global(lambda) only trivially at controls.

    Expected pattern:
      lambda = 0  : intersection-dim = 1 (the global constant — partial)
      lambda = 4,8,10 : intersection-dim = 0 (empty)
      lambda = 12 : intersection-dim = 10 (full)
    """
    by_lambda = {c["lambda"]: c for c in result.negative_controls}
    assert by_lambda[0.0]["intersection_dim"] == 1
    assert by_lambda[4.0]["intersection_dim"] == 0
    assert by_lambda[8.0]["intersection_dim"] == 0
    assert by_lambda[10.0]["intersection_dim"] == 0


def test_exact_rational_certificate_passes(result):
    """The sympy ℚ-rational nullspace certificate must report zero off-coset
    residual components — i.e. the lift identity holds *exactly* over ℚ, not
    just to floating-point precision.
    """
    cert = result.exact_certificate
    assert cert.get("certified") is True, f"exact certificate did not pass: {cert}"
    assert cert["total_local_dim"] == 10
    assert cert["nonzero_components_total"] == 0


def test_exact_minus_one_acts_trivially(result):
    """-1 ∈ 2I acts trivially on E_V600(12) exactly over ℚ. This is what
    allows the right-multiplication action of 2I on the λ=12 eigenspace to
    descend through {±1} and factor through A_5.
    """
    cert = result.exact_certificate
    assert cert["a5_exact"]["minus_one_acts_trivially"] is True


def test_exact_2I_invariance_of_lifted_and_residual(result):
    """E_lifted and E_residual are each preserved by right multiplication
    by every element of 2I, verified exactly over ℚ.

    This is the certificate that supports Lemma 3.5 in the note: for every
    h ∈ 2I, the exact ℚ-rational columns of P_h B lie in the column span
    of B, where B is the rational basis of E_lifted (resp. E_residual).
    Without this lemma the A_5 action on each summand is not well-defined.
    """
    a5 = result.exact_certificate["a5_exact"]
    assert a5["lifted_2I_invariant"] is True
    assert a5["residual_2I_invariant"] is True


def test_exact_a5_irrep_multiplicities_are_integers(result):
    """Exact ℚ-rational character computation must give integer multiplicities."""
    mult = result.exact_certificate["a5_exact"]["irrep_multiplicities"]
    assert mult["lifted"]["5"] == 2
    assert mult["residual"]["5"] == 3
    assert mult["global"]["5"] == 5
    # all other irreps zero
    for name in ("lifted", "residual", "global"):
        for irrep in ("1", "3", "3p", "4"):
            assert mult[name][irrep] == 0


def test_exact_integer_spectrum_multiplicities(result):
    """Every integer eigenvalue of L_C and L_V600 (plus λ=4,8,10 for L_V600
    even though they are not eigenvalues — exactly certifying dim 0 there)
    has its multiplicity certified by exact ℚ-rational sympy nullspace
    computation, matching the multiplicities listed in
    spectral_bridge_note.tex and the rows of the selectivity table.
    """
    isc = result.exact_certificate["integer_spectrum_check"]
    assert isc["local_integer_mults"] == {0: 1, 4: 4, 8: 9, 10: 8, 12: 2}
    assert isc["global_integer_mults"] == {
        0: 1, 4: 0, 8: 0, 9: 16, 10: 0, 12: 25, 14: 36, 15: 16,
    }


def test_exact_a5_class_characters_are_integers(result):
    """The full class-character table must consist of exact integers."""
    cls_chars = result.exact_certificate["a5_exact"]["class_characters"]
    expected = {
        "1A": {"lifted": 10, "residual": 15, "global": 25},
        "2A": {"lifted": 2,  "residual": 3,  "global": 5},
        "3A": {"lifted": -2, "residual": -3, "global": -5},
        "5A": {"lifted": 0,  "residual": 0,  "global": 0},
        "5B": {"lifted": 0,  "residual": 0,  "global": 0},
    }
    for cls, vals in expected.items():
        for name, expected_chi in vals.items():
            assert cls_chars[cls][name] == expected_chi, (
                f"class {cls}, basis {name}: got {cls_chars[cls][name]}, "
                f"expected {expected_chi}"
            )


# --- 6. exact dimensions 10 / 15 / 25 ----------------------------------

def test_exact_lifted_dimension_is_10(result):
    assert result.lambda12_lifted_dim == 10


def test_exact_residual_dimension_is_15(result):
    assert result.lambda12_residual_dim == 15


def test_exact_global_eigenspace_dimension_is_25(result):
    assert result.lambda12_global_dim == 25


# --- 7. A_5 representation strings present -----------------------------

def test_lifted_decomposition_is_two_Y5(result):
    assert result.summary_strings["E_lifted"] == "2*Y5"


def test_residual_decomposition_is_three_Y5(result):
    assert result.summary_strings["E_residual"] == "3*Y5"


def test_global_decomposition_is_five_Y5(result):
    assert result.summary_strings["E_global"] == "5*Y5"


# --- 8. README exists and contains no overclaiming terms ---------------

def test_readme_exists_and_no_overclaiming():
    readme = README_PATH
    assert readme.exists(), f"missing: {readme}"
    text = readme.read_text().lower()
    forbidden = [
        "proves consciousness",
        "derives the standard model",
        "proves cosmology",
        "proves the universe is a hypersphere",
        "undeniable proof",
        "impossible to refute",
        "final theory",
    ]
    for phrase in forbidden:
        assert phrase not in text, f"README contains overclaim: '{phrase}'"


# --- 9. Construction integrity -----------------------------------------

def test_120_vertices_and_5x24_partition(result):
    assert result.n_vertices == 120
    assert result.coset_sizes == [24] * 5


def test_global_graph_is_12_regular_with_720_edges(result):
    assert result.global_degree == 12
    assert result.global_edges == 720


def test_local_graphs_are_8_regular_with_96_edges(result):
    assert result.local_degrees == [8] * 5
    assert result.local_edges == [96] * 5


# --- 10. Other artifacts ------------------------------------------------

def test_repro_log_present():
    p = OUTPUTS / "wo008_reproducibility_log.txt"
    assert p.exists()
    text = p.read_text()
    assert "WO-008" in text


def test_spectral_bridge_note_md_present():
    p = DOCS / "spectral_bridge_note.md"
    assert p.exists(), f"missing: {p}"
    text = p.read_text()
    assert "lambda" in text.lower() or "λ" in text
    assert "24" in text and "600" in text


def test_spectral_bridge_note_tex_present():
    p = DOCS / "spectral_bridge_note.tex"
    assert p.exists(), f"missing: {p}"
    text = p.read_text()
    assert "\\documentclass" in text
    assert "24" in text and "600" in text


# Marketing copy (public_post_draft.md, webpage_copy.md) is an internal
# artifact of the main VFD research repo, not a load-bearing part of the
# Note.  When present (main repo) we sanity-check it; in the public bundle
# the files are absent and the checks are skipped.
@pytest.mark.skipif(not (DOCS / "public_post_draft.md").exists(),
                    reason="public_post_draft.md not present in this layout")
def test_public_post_draft_present_if_main_repo():
    assert (DOCS / "public_post_draft.md").exists()


@pytest.mark.skipif(not (DOCS / "webpage_copy.md").exists(),
                    reason="webpage_copy.md not present in this layout")
def test_webpage_copy_present_if_main_repo():
    assert (DOCS / "webpage_copy.md").exists()


def test_figures_present():
    figdir = DOCS / "figures"
    for name in (
        "bridge_architecture.svg",
        "coset_decomposition.svg",
        "eigenspace_decomposition.svg",
        "negative_controls.svg",
    ):
        assert (figdir / name).exists(), f"missing figure: {figdir / name}"
