# Changelog

All notable changes to this programme are recorded here. The
programme is a series of reproducible exact-arithmetic notes on
the finite geometry and spectral structure of the 600-cell; each
formal claim is certified by exact rational arithmetic (ℚ or
ℚ(√5)) in the accompanying code.

This repository follows [Semantic Versioning](https://semver.org/)
at the bundle level — MAJOR for new published notes, MINOR for
substantial revisions of existing notes, PATCH for typo /
metadata fixes.

---

## [2.1.0] — 2026-05-13

### Added

- **Paper 3 — From Schläfli Decomposition to Spectral Bridge: A
  First Explicit Closure-Projection Channel in V₆₀₀.**  A
  five-page synthesis / programme map that connects Papers 1 and 2
  into a single picture and defines *closure-projection channel*
  cautiously as the local-to-global compatibility relation the
  two prior papers together exhibit at λ=12.

  Files: `docs/03-closure-projection-channel/v600_first_closure_projection_channel.{tex,pdf,md}`,
  four supporting figures, public-post draft, and long-form web copy.

  Paper 3 introduces **no new mathematical claim** and **no new
  computational certificate**; it depends entirely on the existing
  exact certificates in Papers 1 and 2.

### Changed

- Top-level `README.md` now describes the programme as a
  three-note narrative.
- `CITATION.cff` bumped to 2.1.0; Paper 3 added as a referenced
  report alongside Papers 1 and 2.

## [2.0.0] — 2026-05-13

### Added

- **Paper 1 — The Schläfli Decomposition of the 600-cell: Five
  24-cell Cosets and the Induced A₅ Action.**  Establishes V₆₀₀ =
  2I, V₂₄ = 2T, the five-coset partition, the intrinsic 24-cell
  structure on each coset, and the induced A₅ coset action.  Five
  exact ℚ(√5) certificates, 23 tests.

  Files: `docs/01-schlafli-decomposition/`,
  `closure_transform_engine/decomposition.py`,
  `closure_transform_engine/examples/run_wo007_schlafli_decomposition.py`,
  `closure_transform_engine/tests/test_wo007_schlafli_decomposition.py`.

### Changed

- The repository, originally a single-paper artifact for the
  spectral bridge, is now a two-paper programme.  All Paper 2
  files moved into `docs/02-spectral-bridge/` via `git mv` (history
  preserved).
- Top-level `README.md` rewritten as an umbrella narrative.
- Runners and tests are now path-aware across three layouts (main
  research repo, combined public bundle, legacy single-paper
  bundle).

## [1.0.1] — 2026-05-12

### Added

- Author and contact information added to the Note
  (`docs/spectral_bridge_note.{tex,pdf,md}`).

## [1.0.0] — 2026-05-12

### Added

- Initial release: **The 24--600 Spectral Bridge: A selective λ=12
  eigenspace embedding from five 24-cell cosets into the 600-cell.**
  λ=12 lift identity, integer-eigenvalue multiplicities,
  {±1}-triviality, 2I-invariance of E_lifted / E_residual, exact
  integer A₅ class characters, irrep decomposition
  `2·Y₅ + 3·Y₅ = 5·Y₅`.  Hostile-reviewed over 9 codex rounds
  before public release.
