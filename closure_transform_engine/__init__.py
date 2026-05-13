"""closure_transform_engine: keystone artifact for the 24-600 spectral bridge.

This package wraps a single reproducible geometric result around the
existing vfd_v600 icosian/quaternion infrastructure: the lambda = 12
Laplacian eigenspace of each of the five 24-cell cosets in the 600-cell
lifts by zero-extension into the lambda = 12 eigenspace of the full
600-cell graph.

The result is narrow by design. It does not depend on any physics
interpretation; it is a geometric/spectral fact about the 600-cell.

See `docs/keystone_24_600/` for the public artifact and
`closure_transform_engine.keystone` for the computational API.
"""

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
