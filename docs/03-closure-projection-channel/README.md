# Paper 3 — From Schläfli Decomposition to Spectral Bridge

*A First Explicit Closure-Projection Channel in V₆₀₀.*

> **Status:** Synthesis note / programme map. Not peer-reviewed.

This note introduces no new load-bearing mathematical claim. It
synthesises Papers 1 and 2 into a single picture and defines the
term *closure-projection channel* modestly. It carries no formal
certificate of its own; all formal claims relied upon are certified
in Papers 1 and 2.

## Files

- 📄 **Paper PDF:** [`v600_first_closure_projection_channel.pdf`](v600_first_closure_projection_channel.pdf)
- 📝 **Markdown rendering:** [`v600_first_closure_projection_channel.md`](v600_first_closure_projection_channel.md)
- 🔧 **LaTeX source:** [`v600_first_closure_projection_channel.tex`](v600_first_closure_projection_channel.tex)
- 🖼️ **Figures** (4 SVGs): [`figures/`](figures/)
  - `the_ladder.svg` — full pipeline (2I → 5 cosets → local Laplacians → λ=12 lift → global sector)
  - `local_vs_global_shell.svg` — the two d² shells side by side
  - `sector_to_global.svg` — five local sectors feeding one global sector
  - `a5_rep_split.svg` — `5·Y₅ = 2·Y₅ ⊕ 3·Y₅`
- 📰 **Public copy:**
  - [`public_post_draft.md`](public_post_draft.md) — short X/Mastodon thread draft
  - [`webpage_copy.md`](webpage_copy.md) — long-form site copy

## One-paragraph summary

The Schläfli decomposition (Paper 1) supplies *local sectors*: five
disjoint 24-cell cosets inside V₆₀₀. The 24–600 spectral bridge
(Paper 2) supplies an *exact compatibility channel*: the local λ=12
eigenspaces zero-extend identically into the global λ=12 eigenspace.
Read together, they exhibit a pattern — a local shell structure
inside the 600-cell does not merely sit inside the global geometry;
at one specific spectral value it couples *exactly* to a global
invariant sector. The note names this pattern a *closure-projection
channel* (Definition 1) and explains carefully that the term is
interpretive: the formal content is the combination of the two
exact certificates in Papers 1 and 2, not a new theorem.

## Why the bridge is non-trivial

Papers 1 and 2 use **different Laplacians on different graphs**: the
local Laplacian acts on the 24-cell edge graph at d²=1, and the
global Laplacian acts on the 600-cell edge graph at d²=(3−√5)/2.
These graphs share *no* edges. The lift identity is therefore not a
restriction theorem — it is a **cancellation theorem**: a
finite-rank linear condition that turns out to be satisfied
identically on the local λ=12 eigenspace, and nowhere else in the
local spectrum.

## Where it sits in the programme

This is **Paper 3** — the synthesis / reader's map.

- ← **Paper 1** ([`../01-schlafli-decomposition/`](../01-schlafli-decomposition/))
  — the foundation.
- ← **Paper 2** ([`../02-spectral-bridge/`](../02-spectral-bridge/))
  — the spectral consequence.

The note ends by stating the question the next paper will attempt
to formalise: *given a local sector, a local operator, a global
operator, and a projection/lift map, when does a local eigenspace
contribute to a global invariant sector?*
