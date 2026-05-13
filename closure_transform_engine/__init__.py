"""closure_transform_engine: VFD finite-geometry / spectral artifacts.

Two reproducible artifacts live in this package:

  * ``decomposition`` — the Schläfli decomposition of the 600-cell into
    five right cosets of the binary tetrahedral subgroup 2T, with the
    induced A_5 action.  See ``docs/schlafli_24_600/`` and the public
    bundle at ``release-bundles/the-600-cell-five-24-cell-cosets/``.

  * ``keystone`` — the 24-600 spectral bridge: the selective λ=12
    eigenspace lift from the five 24-cell cosets into the full 600-cell
    Laplacian.  See ``docs/keystone_24_600/`` and the public bundle at
    ``release-bundles/the-24-600-spectral-bridge/``.

Both rely only on the vendored ``vfd_v600`` icosian package for exact
Q(sqrt 5) arithmetic, and produce ℚ-rational (or ℚ(√5)-rational)
certificates for every formal claim.
"""

# Foundation paper — Schläfli decomposition.
from .decomposition import (
    DecompositionResult,
    right_cosets,
    distance_shells_v600,
    certify_v600_construction,
    certify_v24_subgroup,
    certify_right_cosets,
    certify_each_coset_is_24cell,
    certify_a5_coset_action,
    run_decomposition,
)

# Spectral bridge paper.
from .keystone import (
    KeystoneResult,
    coset_partition,
    build_d1_graphs,
    lambda12_lift,
    a5_character_decomposition,
    exact_certify_lambda12_lift,
    exact_a5_characters,
    exact_integer_spectrum_check,
    run_keystone,
)

__all__ = [
    # decomposition
    "DecompositionResult",
    "right_cosets",
    "distance_shells_v600",
    "certify_v600_construction",
    "certify_v24_subgroup",
    "certify_right_cosets",
    "certify_each_coset_is_24cell",
    "certify_a5_coset_action",
    "run_decomposition",
    # keystone / spectral bridge
    "KeystoneResult",
    "coset_partition",
    "build_d1_graphs",
    "lambda12_lift",
    "a5_character_decomposition",
    "exact_certify_lambda12_lift",
    "exact_a5_characters",
    "exact_integer_spectrum_check",
    "run_keystone",
]
