# The 24--600 Spectral Bridge

**A three-note reproducible exact-arithmetic programme on the finite
geometry and spectral structure of the 600-cell.**

> **Status:** Reproducible mathematical notes / exact computational
> certificates (Papers 1–2) and a synthesis / programme map (Paper 3).
> Not peer-reviewed.

The repository hosts three notes meant to be read in order. The
first establishes the finite-geometric decomposition of the 600-cell
into five 24-cells; the second proves a spectral consequence of that
decomposition; the third explains their combined role as a first
local-to-global closure-projection channel.

---

## Read in order

### Paper 1 — The foundation

**The Schläfli Decomposition of the 600-cell: Five 24-cell Cosets
and the Induced A₅ Action**

📄 [`docs/01-schlafli-decomposition/schlafli_decomposition.pdf`](docs/01-schlafli-decomposition/schlafli_decomposition.pdf) ·
[markdown](docs/01-schlafli-decomposition/schlafli_decomposition.md) ·
[LaTeX source](docs/01-schlafli-decomposition/schlafli_decomposition.tex)

The 600-cell vertex set is the binary icosahedral group 2I (120 unit
icosian quaternions). The binary tetrahedral group 2T is the σ-fixed
subgroup of 2I, of index 5; its five right cosets partition V₆₀₀
into five disjoint 24-element subsets, each carrying the intrinsic
24-cell distance structure. The induced action of 2I on these five
cosets factors through {±1} to a faithful transitive A₅. **Every
claim is certified by exact ℚ(√5) arithmetic; no floating-point
computation enters any formal claim.**

### Paper 2 — The spectral consequence

**The 24--600 Spectral Bridge: A selective λ=12 eigenspace embedding
from five 24-cell cosets into the 600-cell**

📄 [`docs/02-spectral-bridge/spectral_bridge_note.pdf`](docs/02-spectral-bridge/spectral_bridge_note.pdf) ·
[markdown](docs/02-spectral-bridge/spectral_bridge_note.md) ·
[LaTeX source](docs/02-spectral-bridge/spectral_bridge_note.tex)

Building on Paper 1, the 2-dimensional λ=12 eigenspace of each local
24-cell Laplacian zero-extends *exactly* into the 25-dimensional
λ=12 eigenspace of the full 600-cell Laplacian. The lift is
*selective*: no other local Laplacian eigenvalue admits a
non-trivial lift. Under the induced A₅ action on the five cosets,
the global λ=12 eigenspace decomposes as `5·Y₅ = 2·Y₅ ⊕ 3·Y₅`, with
exact integer characters. **Every claim is certified by exact
ℚ-rational (or ℚ(√5)-rational) arithmetic.**

### Paper 3 — The synthesis

**From Schläfli Decomposition to Spectral Bridge: A First
Closure-Projection Channel in V₆₀₀**

📄 [`docs/03-closure-projection-channel/v600_first_closure_projection_channel.pdf`](docs/03-closure-projection-channel/v600_first_closure_projection_channel.pdf) ·
[markdown](docs/03-closure-projection-channel/v600_first_closure_projection_channel.md) ·
[LaTeX source](docs/03-closure-projection-channel/v600_first_closure_projection_channel.tex)

A short (5-page) reader's map that synthesises Papers 1 and 2 into a
single picture: a local 24-cell shell structure couples *exactly* to
a global 600-cell invariant sector at λ=12, with A₅-stable
decomposition `2·Y₅ ⊕ 3·Y₅ = 5·Y₅`. Defines "closure-projection
channel" cautiously as a name for the local-to-global compatibility
relation the two prior papers together exhibit — not as a new
theorem. The note introduces no new mathematical claim and no
physics claim; it carries no formal certificate of its own and
depends entirely on the exact certificates in Papers 1 and 2.

---

## Reproduce

```bash
pip install -r requirements.txt          # or: pip install -e .

# Paper 1 — foundation (under a minute):
python closure_transform_engine/examples/run_wo007_schlafli_decomposition.py

# Paper 2 — spectral bridge (~1.5 minutes):
python closure_transform_engine/examples/run_wo008_keystone.py

# Re-assert every certificate-level claim as a test (~2.5 minutes):
pytest closure_transform_engine/tests/
```

Both runners write frozen output artifacts into the corresponding
`docs/0?-*/outputs/` directory, including a reproducibility log that
records UTC timestamp, environment, and commit hash of the source
revision used.

Requires Python ≥ 3.8 with NumPy, SymPy, SciPy (the `pyproject.toml`
specifies minimum versions; developed on Python 3.8.10, NumPy 1.24.4,
SymPy 1.13.3, SciPy 1.10.1).

---

## What each runner verifies (one-line summaries)

### Paper 1 — five certificates

| § | Theorem | Certificate routine |
|---|---|---|
| 2 | V₆₀₀ is a group of order 120 (icosian model of 2I) | `certify_v600_construction` |
| 3 | V₂₄ = σ-fixed subset is a subgroup of index 5 (Hurwitz-unit model of 2T) | `certify_v24_subgroup` |
| 4 | Five right cosets of 2T partition V₆₀₀, each size 24 | `certify_right_cosets` |
| 5 | Each coset is isometric to the 24-cell; no intra-coset V₆₀₀-edge | `certify_each_coset_is_24cell` |
| 6 | Induced 2I-action on five cosets has kernel {±1} and image A₅ | `certify_a5_coset_action` |

### Paper 2 — five certificates (plus the floating-point cross-check)

| § | Claim | Certificate routine |
|---|---|---|
| 3 | λ=12 lift identity (Theorem 1) | `exact_certify_lambda12_lift` |
| 2/4 | Integer-eigenvalue multiplicities of L_C and L_V600 (incl. λ=4,8,10 nullity = 0) | `exact_integer_spectrum_check` |
| 5 | Lemma 3 ({±1}-triviality on E_V600(12)) | `exact_a5_characters.minus_one_acts_trivially` |
| 5 | Lemma 3.5 (2I-invariance of E_lifted and E_residual) | `exact_a5_characters.{lifted,residual}_2I_invariant` |
| 5 | Theorem 4 + Corollary 5 (exact integer A₅ characters and irrep decomposition) | `exact_a5_characters` |
| 3 | Floating-point cross-check `‖L₆₀₀ v − 12 v‖ ≤ 6.2e-15` | `lambda12_lift` (NumPy) |

The NumPy floating-point computation is used **only** for the
independent cross-check; it does not enter any formal claim.

---

## Repository layout

```
README.md                        this file (umbrella narrative)
LICENSE                          MIT
CITATION.cff                     citation metadata for all three notes
requirements.txt, pyproject.toml  dependency pins (Python ≥ 3.8)
.gitignore

closure_transform_engine/        the computation (Papers 1 + 2)
  decomposition.py               Paper 1 certificates
  keystone.py                    Paper 2 certificates
  examples/
    run_wo007_schlafli_decomposition.py   Paper 1 runner
    run_wo008_keystone.py                  Paper 2 runner
  tests/
    test_wo007_schlafli_decomposition.py  Paper 1 tests (23)
    test_wo008_keystone_artifact.py       Paper 2 tests (28)

vfd_v600/                        vendored exact Q(sqrt 5) icosian-quaternion
                                 package (the only mathematical dependency)

docs/
  01-schlafli-decomposition/     Paper 1 artifacts
    schlafli_decomposition.{tex,pdf,md}
    figures/    coset_partition, group_inclusion, distance_shells, coset_action
    outputs/    wo007_summary.json, cosets.json, coset_action_table.csv,
                distance_shells.csv, reproducibility_log.txt
  02-spectral-bridge/            Paper 2 artifacts
    spectral_bridge_note.{tex,pdf,md}
    figures/    bridge_architecture, coset_decomposition,
                eigenspace_decomposition, negative_controls
    outputs/    wo008_summary.json, verdict_table.csv,
                reproducibility_log.txt
  03-closure-projection-channel/  Paper 3 — synthesis / programme map
    v600_first_closure_projection_channel.{tex,pdf,md}
    figures/    the_ladder, local_vs_global_shell, sector_to_global, a5_rep_split
    public_post_draft.md, webpage_copy.md
```

---

## What this programme claims and does not claim

It claims:

- A reproducible exact-arithmetic certificate that V₆₀₀ partitions
  into five isometric 24-cells (the Schläfli decomposition);
- A reproducible exact-arithmetic certificate that under the induced
  2I-action on these five cosets, the kernel is {±1} and the image
  is A₅;
- A reproducible exact-arithmetic certificate that the local 24-cell
  λ=12 eigenspaces zero-extend exactly into the global 600-cell λ=12
  eigenspace, and that no other local Laplacian eigenvalue admits
  a non-trivial lift;
- An exact A₅-equivariant irrep decomposition
  `E_V600(12) = 2·Y₅ ⊕ 3·Y₅`.

It does **not** claim to derive the Standard Model, particle masses,
cosmology, black holes, the Riemann hypothesis, or consciousness,
and it does not subsume or replace any other framework.

---

## Provenance

The two papers stand on a single common foundation:

- `vfd_v600/` (vendored, unchanged) — exact ℚ(√5) icosian quaternion
  arithmetic, the 120-vertex V₆₀₀ vertex set, the σ-Galois twist,
  and the σ-fixed 24-cell V₂₄.
- The classical fact that 2T has index 5 in 2I and that the cosets
  form five inscribed 24-cells is Coxeter / Conway–Sloane; the
  precise references are in each paper's bibliography.

Each paper is self-contained for verification: no external
repository is required to check any formal claim.

## Citation

See `CITATION.cff` for the umbrella citation; each paper's
bibliography also gives detailed references for its own dependencies.

## Changelog

Publication history is in [`CHANGELOG.md`](CHANGELOG.md).

## License

MIT — see `LICENSE`.
