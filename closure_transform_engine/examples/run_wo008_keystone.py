"""WO-008 keystone runner: 24-600 lambda=12 spectral bridge.

Single entry point. Prints the public summary, writes JSON / CSV / log
artifacts under the docs outputs directory, and exits with status 0
on a positive verdict.

Run:
    python closure_transform_engine/examples/run_wo008_keystone.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
from pathlib import Path

# Allow running as a script: ensure repo root is on sys.path
_HERE = Path(__file__).resolve()
_REPO = _HERE.parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from closure_transform_engine.keystone import run_keystone, format_decomposition


def _resolve_out_dir() -> Path:
    """Locate the docs outputs directory in either repo layout.

    Main VFD research repo: docs/keystone_24_600/outputs/.
    Standalone reproduction bundle: docs/outputs/.
    """
    candidates = [
        _REPO / "docs" / "keystone_24_600" / "outputs",
        _REPO / "docs" / "outputs",
    ]
    for c in candidates:
        if c.parent.exists():  # the docs dir exists
            return c
    return candidates[0]


OUT_DIR = _resolve_out_dir()


def emit(line: str, log_lines):
    print(line)
    log_lines.append(line)


def main(argv=None) -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log_lines = []

    t0 = time.time()
    emit("=== WO-008 Keystone: 24-600 lambda=12 Spectral Bridge ===", log_lines)
    emit("", log_lines)
    emit("Constructing 600-cell...", log_lines)
    emit("Computing icosian/quaternion 5-coset partition...", log_lines)
    emit("Verifying each coset carries d=1 24-cell structure...", log_lines)
    emit("Computing local 24-cell eigenspaces...", log_lines)
    emit("Computing global 600-cell lambda=12 eigenspace...", log_lines)
    emit("Testing zero-extension lifts...", log_lines)
    emit("", log_lines)

    result = run_keystone()

    emit("Result:", log_lines)
    emit(f"  lambda=12 lift residual: {result.lambda12_max_residual:.3e}", log_lines)
    emit(f"  lifted dimension: {result.lambda12_lifted_dim}", log_lines)
    emit(f"  residual dimension: {result.lambda12_residual_dim}", log_lines)
    emit(f"  global E600(12) dimension: {result.lambda12_global_dim}", log_lines)
    emit("", log_lines)

    emit("Lift behaviour by local lambda (intersection-dim = dim of lift image"
         "  cap  E_global(lambda)):", log_lines)
    for c in result.negative_controls:
        status = (
            "FULL"  if c["intersection_dim"] == c["lifted_dim"] and c["lifted_dim"] > 0 else
            "PARTIAL" if c["intersection_dim"] > 0 else
            "EMPTY"
        )
        emit(
            f"  lambda={c['lambda']:>4}: lifted-dim={c['lifted_dim']:>2}, "
            f"global-dim={c['global_dim']:>2}, intersection-dim={c['intersection_dim']:>2}, "
            f"max-residual={c['max_residual']:.3e}  {status}",
            log_lines,
        )
    emit(
        f"  lambda=  12: lifted-dim={result.lambda12_lifted_dim:>2}, "
        f"global-dim={result.lambda12_global_dim:>2}, "
        f"intersection-dim={result.lambda12_lifted_dim:>2}, "
        f"max-residual={result.lambda12_max_residual:.3e}  FULL",
        log_lines,
    )
    emit("", log_lines)

    if result.exact_certificate:
        cert = result.exact_certificate
        emit("Exact ℚ-rational certificate (sympy, integer Laplacians):", log_lines)
        emit(f"  total local-λ=12 kernel dim: {cert['total_local_dim']}", log_lines)
        emit(f"  nonzero off-coset components after lift: {cert['nonzero_components_total']}", log_lines)
        emit(
            "  lift identity certified: "
            + ("YES — exactly over ℚ" if cert.get("certified") else "NO"),
            log_lines,
        )
        a5 = cert.get("a5_exact", {})
        if a5:
            emit(f"  -1 ∈ 2I acts trivially on E_V600(12): "
                 + ("YES (exact)" if a5["minus_one_acts_trivially"] else "NO"),
                 log_lines)
            emit("  exact A_5 class characters:", log_lines)
            for cls, ch in a5["class_characters"].items():
                emit(f"    {cls}: lifted={ch['lifted']:>4}  residual={ch['residual']:>4}  global={ch['global']:>4}",
                     log_lines)
            emit("  exact A_5 irrep multiplicities:", log_lines)
            for name in ("lifted", "residual", "global"):
                terms = []
                for irrep, m in a5["irrep_multiplicities"][name].items():
                    if m != 0:
                        terms.append(f"{m}·{irrep}")
                emit(f"    {name}: {' + '.join(terms) if terms else '0'}", log_lines)
        emit("", log_lines)

    emit("Representation (A_5 = 2I/{+/-1} acting on the 5 cosets):", log_lines)
    emit(f"  E_lifted   ~= {result.summary_strings['E_lifted']}", log_lines)
    emit(f"  E_residual ~= {result.summary_strings['E_residual']}", log_lines)
    emit(f"  E600(12)  ~= {result.summary_strings['E_global']}", log_lines)
    emit("", log_lines)

    emit(f"Verdict:", log_lines)
    emit(f"  {result.verdict}", log_lines)
    emit("", log_lines)

    # Write JSON summary
    summary_path = OUT_DIR / "wo008_summary.json"
    with summary_path.open("w") as f:
        json.dump(result.to_dict(), f, indent=2)
    emit(f"Wrote: {summary_path.relative_to(_REPO)}", log_lines)

    # Write verdict CSV
    csv_path = OUT_DIR / "wo008_verdict_table.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "lambda",
            "lift_image_dim",
            "global_eigenspace_dim",
            "intersection_dim",
            "max_lift_residual",
            "verdict",
        ])
        rows = list(result.negative_controls)
        rows.append({
            "lambda": 12.0,
            "lifted_dim": result.lambda12_lifted_dim,
            "global_dim": result.lambda12_global_dim,
            "intersection_dim": result.lambda12_lifted_dim,
            "max_residual": result.lambda12_max_residual,
        })
        for c in rows:
            if c["intersection_dim"] == c["lifted_dim"] and c["lifted_dim"] > 0:
                verdict = "FULL"
            elif c["intersection_dim"] > 0:
                verdict = "PARTIAL"
            else:
                verdict = "EMPTY"
            w.writerow([
                c["lambda"],
                c["lifted_dim"],
                c["global_dim"],
                c["intersection_dim"],
                f"{c['max_residual']:.6e}",
                verdict,
            ])
    emit(f"Wrote: {csv_path.relative_to(_REPO)}", log_lines)

    # Write reproducibility log
    log_path = OUT_DIR / "wo008_reproducibility_log.txt"
    # Record UTC timestamp + commit hash (best-effort; works in any git repo).
    import datetime
    import subprocess
    run_utc = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    try:
        commit_hash = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(_REPO),
            capture_output=True, text=True, timeout=5,
        ).stdout.strip() or "(no git repo / no commits)"
    except Exception:
        commit_hash = "(git unavailable)"
    try:
        commit_dirty = subprocess.run(
            ["git", "status", "--porcelain"], cwd=str(_REPO),
            capture_output=True, text=True, timeout=5,
        ).stdout.strip()
    except Exception:
        commit_dirty = ""
    dirty_flag = " (uncommitted changes present)" if commit_dirty else ""

    with log_path.open("w") as f:
        f.write(f"# WO-008 keystone run log\n")
        f.write(f"# run_utc:  {run_utc}\n")
        f.write(f"# elapsed:  {time.time() - t0:.2f} s\n")
        f.write(f"# python:   {sys.version.splitlines()[0]}\n")
        f.write(f"# repo:     {_REPO}\n")
        f.write(f"# commit:   {commit_hash}{dirty_flag}\n\n")
        for line in log_lines:
            f.write(line + "\n")
    print(f"Wrote: {log_path.relative_to(_REPO)}")

    return 0 if result.verdict == "SELECTIVE_SPECTRAL_BRIDGE" else 1


if __name__ == "__main__":
    sys.exit(main())
