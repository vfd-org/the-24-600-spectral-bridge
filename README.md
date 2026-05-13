# The 24--600 Spectral Bridge

**A reproducible exact-arithmetic spectral-graph result: the selective
λ=12 eigenspace bridge between five 24-cell cosets and the 600-cell.**

> **Status:** Reproducible mathematical note / exact computational
> certificate. Not peer-reviewed.

---

## The result in one paragraph

The 600-cell (the regular 4-polytope with 120 vertices) has a natural
five-fold structure: its vertex set is the binary icosahedral group
2I, the binary tetrahedral group 2T sits inside it with index 5, and
the right cosets of 2T partition the 600-cell into **five disjoint
24-cells**. On each 24-cell coset, the nearest-neighbour ("d=1")
graph is the 24-cell edge graph (8-regular, 96 edges); on the full
vertex set, the d=1 graph is the 600-cell edge graph (12-regular,
720 edges). The two graphs share **no edges** — the local 24-cell
edge length is the next distance shell up. We prove that the
2-dimensional λ=12 eigenspace of each local 24-cell Laplacian
zero-extends into the 25-dimensional λ=12 eigenspace of the full
600-cell Laplacian, exactly over ℚ, and that **no other local
Laplacian eigenvalue admits a non-trivial lift**. Under the induced
A₅ = 2I/{±1} action on the five cosets, the 25-dim global λ=12
eigenspace decomposes as `2·Y₅ ⊕ 3·Y₅` where Y₅ is the 5-dim
irreducible representation of A₅, with exact integer characters.

## Why the comparison is non-trivial

A function supported on a single coset is generically *not* a
`L₆₀₀`-eigenvector: applying `L₆₀₀` couples the coset to vertices
outside it (via 600-cell edges that cross between cosets). The lift
identity says that for the local λ=12 eigenspace this off-coset
coupling cancels *identically* — a finite-rank linear condition that
turns out to be satisfied. The selectivity table says no other local
eigenspace satisfies it. That is what makes λ=12 a genuine
local-to-global spectral channel, not a trivial
restriction-and-extension.

## The formal claims (all certified exactly over ℚ)

The proof package is a set of exact rational-arithmetic certificates,
not floating-point evidence. Every one of the following is checked in
ℚ (or ℚ(√5) where character values are irrational) by
`closure_transform_engine/keystone.py` and re-asserted by the test
suite:

| Claim | Certificate |
|---|---|
| **Theorem 1** — λ=12 lift identity: `liftᵢ(E_{Cᵢ}(12)) ⊂ E₆₀₀(12)` for every coset | `exact_certify_lambda12_lift` — integer Laplacians, ℚ-rational nullspace bases, exact zero off-coset residuals in all 10 lifted basis vectors |
| **Selectivity table** — `dim E₆₀₀(λ)` for λ ∈ {0, 4, 8, 9, 10, 12, 14, 15} | `exact_integer_spectrum_check` — ℚ-rational sympy nullspace at each integer λ |
| **Lemma 3** — `−1 ∈ 2I` acts trivially on `E₆₀₀(12)` | `exact_a5_characters.minus_one_acts_trivially` — componentwise check of `P₋₁ v − v = 0` over ℚ for all 25 rational basis vectors |
| **Lemma 3.5** — `E_lifted` and `E_residual` are 2I-invariant | `exact_a5_characters.{lifted,residual}_2I_invariant` — for every h ∈ 2I and B ∈ {B_lifted, B_residual}, exact rational check of `P_h B − B(BᵀB)⁻¹Bᵀ(P_h B) = 0` |
| **Theorem 4** — exact integer A₅ class characters | exact rational trace, A₅ conjugacy classes enumerated by orbit (verified sizes 1, 15, 20, 12, 12) |
| **Corollary 5** — irrep decomposition `E_lifted ≅ 2·Y₅`, `E_residual ≅ 3·Y₅`, `E₆₀₀(12) ≅ 5·Y₅` | exact ℚ(√5) inner products against the A₅ character table |

Floating-point computation (NumPy) enters only as an *independent
cross-check* in `lambda12_lift` / `run_keystone`, where it agrees
with the exact certificates to ‖`L₆₀₀ v − 12 v`‖ ≤ 6.2 × 10⁻¹⁵.

## Selectivity (full table)

| local λ | lift-image dim | dim E₆₀₀(λ) | intersection dim | verdict |
|---------|---------------:|------------:|-----------------:|---------|
| 0 | 5 | 1 | 1 | PARTIAL — only the global constant (symmetric combination of coset indicators) |
| 4 | 20 | 0 | 0 | E₆₀₀(λ) = {0} (vacuous; λ=4 is not an eigenvalue of `L₆₀₀`) |
| 8 | 45 | 0 | 0 | E₆₀₀(λ) = {0} |
| 10 | 40 | 0 | 0 | E₆₀₀(λ) = {0} |
| **12** | **10** | **25** | **10** | **FULL** |

Among the five local 24-cell Laplacian eigenvalues, **λ=12 is the only
one that gives a full zero-extension** into a corresponding global
600-cell eigenspace. The λ=0 case contributes exactly one direction
(the global constant, reconstructed from the symmetric sum of coset
indicators). The local eigenvalues λ ∈ {4, 8, 10} are not eigenvalues
of `L₆₀₀` at all (certified by exact ℚ-rational nullspace of
`L₆₀₀ − λI`), so the question is vacuous there.

## What this does and does not claim

It claims:
- a reproducible, selective λ=12 spectral channel between the five
  24-cell cosets of the 600-cell and the global 25-dim λ=12 sector;
- an exact direct-sum split of dimensions 10 + 15 = 25;
- an A₅-equivariant decomposition `2·Y₅ + 3·Y₅ = 5·Y₅` with exact
  integer characters.

It does **not** claim to derive the Standard Model, particle masses,
cosmology, black holes, the Riemann hypothesis, or consciousness, and
it does not subsume or replace any other framework. It is a narrow,
inspectable, reproducible geometry result.

## Reproduce it

Requires Python ≥ 3.8 with NumPy, SciPy, SymPy (see
`requirements.txt`; developed on Python 3.8.10, NumPy 1.24.4, SymPy
1.13.3, SciPy 1.10.1). The formal certificates use only exact ℚ
arithmetic via SymPy; NumPy is used only for the floating-point
cross-check.

```bash
pip install -r requirements.txt          # or: pip install -e .

# Run the full protocol (a few minutes): builds the 600-cell from
# exact Q(sqrt 5) icosian quaternions, runs every certificate, and
# writes the frozen output artifacts to docs/outputs/.
python closure_transform_engine/examples/run_wo008_keystone.py

# Re-assert the certificate-level claims as a test (~2 minutes):
pytest closure_transform_engine/tests/test_wo008_keystone_artifact.py
```

Expected console summary:

```text
λ=12 lift residual: ~6.16e-15  (numerical cross-check)
exact ℚ-rational certificate: YES (0 nonzero off-coset components)
-1 ∈ 2I acts trivially on E_V600(12): YES (exact)
E_lifted, E_residual 2I-invariant: YES (exact, all 120 elements)
selectivity:  λ=0 PARTIAL (intersection 1/5)
              λ=4,8,10 EMPTY (intersection 0; dim E_V600 = 0 exactly)
              λ=12 FULL (intersection 10/10)
exact A_5 irrep multiplicities:  E_lifted 2·Y5,  E_residual 3·Y5,  E600(12) 5·Y5
Verdict: SELECTIVE_SPECTRAL_BRIDGE
```

The reproducibility log (`docs/outputs/wo008_reproducibility_log.txt`)
records the run UTC timestamp, elapsed time, Python version,
repository path, and the **commit hash** of the source revision used
(plus a marker if the working tree was dirty), so the frozen
artifacts are anchored to a specific source state.

## Repository layout

```
closure_transform_engine/        the computation
  keystone.py                    construction, lift, exact certificates, A_5 characters
  examples/run_wo008_keystone.py one-command runner
  tests/test_wo008_keystone_artifact.py   the test suite (28 tests; asserts the certificate-level claims)
vfd_v600/                        vendored exact Q(sqrt 5) icosian-quaternion package
                                 (V_600 = 2I, the sigma-Galois twist, the 24-cell V_24,
                                 and the V_24 left action used for the five-coset partition)
docs/
  spectral_bridge_note.{tex,pdf,md}   the technical note (canonical = .tex)
  figures/                       SVG figures
  outputs/                       frozen run artifacts (JSON / CSV / log)
requirements.txt, pyproject.toml, LICENSE, CITATION.cff
```

## Provenance and dependencies

The bundle is self-contained: **the claims of this note do not
require any external repository.** A reader can verify every theorem
by running the two commands above.

For context, the underlying construction connects to a broader
research programme:

- `vfd_v600/` (vendored unchanged in this repo) provides exact ℚ(√5)
  icosian quaternion arithmetic, the 120-vertex V₆₀₀ vertex set, the
  σ-Galois twist, the σ-fixed 24-cell V₂₄, and the V₂₄ left action
  used to obtain the five-coset partition.
- The classical fact that 2T has index 5 in 2I, and that the cosets
  form five inscribed 24-cells, is Coxeter / Conway–Sloane (cited in
  the note's bibliography). A machine-verified treatment within the
  upstream research programme uses the left-coset convention; this
  repo uses right cosets, related by the involution g ↦ g⁻¹.
- The 600-cell and 24-cell graph spectra are also computed
  independently elsewhere in the research programme; the integer
  multiplicities used in the note are additionally certified by the
  exact ℚ-rational nullspace routines included here.

## Citation

See `CITATION.cff`.

## License

MIT — see `LICENSE`.
