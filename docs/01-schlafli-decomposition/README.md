# Paper 1 — The Schläfli Decomposition of the 600-cell

*Five 24-cell Cosets and the Induced A₅ Action.*

> **Status:** Reproducible mathematical note / exact computational
> certificate. Not peer-reviewed.

## Files

- 📄 **Paper PDF:** [`schlafli_decomposition.pdf`](schlafli_decomposition.pdf)
- 📝 **Markdown rendering:** [`schlafli_decomposition.md`](schlafli_decomposition.md)
- 🔧 **LaTeX source:** [`schlafli_decomposition.tex`](schlafli_decomposition.tex)
- 🖼️ **Figures** (4 SVGs): [`figures/`](figures/)
- 📊 **Frozen outputs:** [`outputs/`](outputs/)
  - `wo007_summary.json` — full structured certificate result
  - `wo007_cosets.json` — the five coset vertex-index lists
  - `wo007_coset_action_table.csv` — induced coset permutation per element of 2I
  - `wo007_distance_shells.csv` — the 8 non-zero V₆₀₀ distance shells
  - `wo007_reproducibility_log.txt` — UTC timestamp, environment, commit hash

## One-paragraph summary

The 600-cell vertex set V₆₀₀ is the binary icosahedral group 2I (120
unit icosian quaternions). The binary tetrahedral group 2T sits
inside 2I as the σ-fixed subgroup under the Galois twist
`σ: √5 ↦ −√5`, with index 5. Its five right cosets partition V₆₀₀
into five disjoint 24-element subsets, each isometric to the
standard 24-cell as a metric subspace of S³. The induced action of
2I on these five cosets is transitive, has kernel exactly {±1}, and
image A₅ = 2I/{±1}. **Every claim is certified by exact ℚ(√5)
arithmetic.**

## Reproduce

From the repository root:

```bash
python ../../closure_transform_engine/examples/run_wo007_schlafli_decomposition.py
pytest  ../../closure_transform_engine/tests/test_wo007_schlafli_decomposition.py
```

23 tests covering: V₆₀₀ construction, V₂₄ subgroup, five-coset
partition, 24-cell structure on each coset, A₅ coset action,
reproducibility artifacts.

## Where it sits in the programme

This is **Paper 1** — the finite-geometric foundation. It supplies
the local sectors used downstream.

- → **Paper 2** ([`../02-spectral-bridge/`](../02-spectral-bridge/))
  — the spectral consequence: local λ=12 eigenspaces lift exactly
  into the global λ=12 sector.
- → **Paper 3** ([`../03-closure-projection-channel/`](../03-closure-projection-channel/))
  — synthesis: how Papers 1+2 combine into a first explicit
  closure-projection channel.
