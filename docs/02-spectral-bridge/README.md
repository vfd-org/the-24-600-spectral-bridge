# Paper 2 — The 24--600 Spectral Bridge

*A selective λ=12 eigenspace embedding from five 24-cell cosets into
the 600-cell.*

> **Status:** Reproducible mathematical note / exact computational
> certificate. Not peer-reviewed.

## Files

- 📄 **Paper PDF:** [`spectral_bridge_note.pdf`](spectral_bridge_note.pdf)
- 📝 **Markdown rendering:** [`spectral_bridge_note.md`](spectral_bridge_note.md)
- 🔧 **LaTeX source:** [`spectral_bridge_note.tex`](spectral_bridge_note.tex)
- 🖼️ **Figures** (4 SVGs): [`figures/`](figures/)
- 📊 **Frozen outputs:** [`outputs/`](outputs/)
  - `wo008_summary.json` — full structured certificate result
  - `wo008_verdict_table.csv` — per-λ verdict (FULL / PARTIAL / EMPTY)
  - `wo008_reproducibility_log.txt` — UTC timestamp, environment, commit hash

## One-paragraph summary

Building on Paper 1's five-coset partition of V₆₀₀ = 2I, the
2-dimensional λ=12 eigenspace of each local 24-cell Laplacian
zero-extends *exactly* into the 25-dimensional λ=12 eigenspace of
the full 600-cell Laplacian. The lift is *selective*: at λ ∈ {4,
8, 10} the global eigenspace is `{0}` (certified exactly over ℚ), so
the question is vacuous; at λ=0 only the global constant survives,
through the symmetric sum of coset indicators. Under the induced A₅
action, `E_V600(12) = 2·Y₅ ⊕ 3·Y₅` with exact integer characters.
**Every claim is certified by exact ℚ-rational (or ℚ(√5)-rational)
arithmetic; floating-point computation is included only as an
independent cross-check.**

## Reproduce

From the repository root:

```bash
python ../../closure_transform_engine/examples/run_wo008_keystone.py
pytest  ../../closure_transform_engine/tests/test_wo008_keystone_artifact.py
```

28 tests covering: λ=12 lift identity, integer-eigenvalue
multiplicities, {±1}-triviality, 2I-invariance of E_lifted /
E_residual, exact integer A₅ class characters, irrep
decomposition, runner integration, and a floating-point cross-check
(`‖L₆₀₀ v − 12v‖ ≤ 6.2 × 10⁻¹⁵`).

## Where it sits in the programme

This is **Paper 2** — the spectral consequence of the foundation
geometry.

- ← **Paper 1** ([`../01-schlafli-decomposition/`](../01-schlafli-decomposition/))
  — the foundation: 5 coset partition + A₅ coset action.
- → **Paper 3** ([`../03-closure-projection-channel/`](../03-closure-projection-channel/))
  — synthesis: how Papers 1+2 combine into a first explicit
  closure-projection channel.
