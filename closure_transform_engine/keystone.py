"""Compute the 24-600 lambda=12 spectral bridge.

Pipeline (each step is exact or numerically dominant):

  1. Construct V_600 = 2I (120 icosian quaternions) via vfd_v600.
  2. Partition V_600 into 5 right cosets of V_24 = 2T.
     - One coset = V_24 itself (sigma-fixed sublattice).
     - Four cosets = the V_24 left-mult orbits on sigma-mobile (96 vertices).
  3. Build the d=1 nearest-neighbour graph on V_600
     (12-regular, 720 edges) and the d=1 graph on each V_24 coset
     (8-regular, 96 edges, isomorphic to the 24-cell edge graph).
  4. Diagonalise the Laplacians.  Extract the local lambda = 12
     eigenspace (2-dim per coset) and global lambda = 12 eigenspace
     (25-dim).
  5. Zero-extend each local eigenvector into R^120 and check
     ||L_600 v - 12 v|| (must be ~ machine epsilon).
  6. Run the same protocol at lambda in {0, 4, 8, 10} as negative
     controls.
  7. Compute the A_5 character of E_lifted (10-dim) and E_residual
     (15-dim = E_600(12) \\ E_lifted), where A_5 = 2I/{+/-1} acts on
     the 5 cosets by right multiplication, and decompose against the
     standard A_5 character table.

All numerical operations use numpy; the underlying vertex and coset
construction is exact (Q(sqrt 5) arithmetic in vfd_v600).
"""
from __future__ import annotations

import math
import sys
import os
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np

# Locate the vfd_v600 icosian package.  It may already be importable
# (e.g. installed, or on PYTHONPATH, or vendored at the repo root next to
# this package); otherwise fall back to known in-tree locations.  This
# keeps the module working both inside the main VFD research repo and
# inside the self-contained reproduction bundle.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
try:  # already importable?
    import vfd_v600  # type: ignore  # noqa: F401
except ImportError:
    for _cand in (
        os.path.join(_REPO),                                            # vendored: <root>/vfd_v600/
        os.path.join(_REPO, "papers", "v600-programme", "lib"),         # in-tree main repo
    ):
        if os.path.isdir(os.path.join(_cand, "vfd_v600")):
            if _cand not in sys.path:
                sys.path.insert(0, _cand)
            break

from vfd_v600.group import build_state  # type: ignore
from vfd_v600.symmetry import mobile_orbit_under_v24  # type: ignore
from vfd_v600.sigma import sigma_fixed_set  # type: ignore
from vfd_v600.icosian import qq_distance_sq, qq_mul, qq_key  # type: ignore


# ----- exact distance utilities -----------------------------------------

_ZERO_D2 = (Fraction(0), Fraction(0))


def _q5_to_float(p) -> float:
    return float(p[0]) + float(p[1]) * math.sqrt(5.0)


# ----- 5-coset partition -------------------------------------------------

def coset_partition(state) -> List[Tuple[int, ...]]:
    """Return the 5 V_24 right-cosets in V_600 as sorted tuples of indices.

    Coset 0 is V_24 (sigma-fixed). Cosets 1..4 are the V_24 left-mult
    orbits on the 96 sigma-mobile vertices. The set-union of all five
    is exactly V_600, and every coset has size 24.
    """
    sf = sigma_fixed_set(state)
    mobile = mobile_orbit_under_v24(state)
    cosets = [tuple(sorted(sf))] + [tuple(sorted(o)) for o in mobile]
    assert len(cosets) == 5
    assert all(len(c) == 24 for c in cosets)
    flat: set = set()
    for c in cosets:
        flat.update(c)
    assert flat == set(range(state["n"]))
    return cosets


# ----- d=1 graph + Laplacians -------------------------------------------

def build_d1_graphs(state, cosets) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Build the V_600 d=1 Laplacian and per-coset 24-cell d=1 Laplacians.

    The V_600 d=1 graph uses the global shortest non-zero distance
    (shell 1, d^2 = (3-sqrt 5)/2). Each per-coset 24-cell d=1 graph
    uses that coset's internal shortest non-zero distance (d^2 = 1
    for every coset, identical because V_24 cosets are isometric).
    """
    verts = state["verts"]
    n = state["n"]

    # Distinct non-zero d^2 values on V_600
    all_d2: set = set()
    for i in range(n):
        for j in range(n):
            d = qq_distance_sq(verts[i], verts[j])
            if d == _ZERO_D2:
                continue
            all_d2.add(d)
    sorted_d2 = sorted(all_d2, key=_q5_to_float)
    d1_global = sorted_d2[0]

    A_global = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if qq_distance_sq(verts[i], verts[j]) == d1_global:
                A_global[i, j] = 1.0
    L_global = np.diag(A_global.sum(axis=1)) - A_global

    L_local: List[np.ndarray] = []
    for c in cosets:
        sub_d2: set = set()
        for i in c:
            for j in c:
                if i == j:
                    continue
                d = qq_distance_sq(verts[i], verts[j])
                if d == _ZERO_D2:
                    continue
                sub_d2.add(d)
        d1_local = sorted(sub_d2, key=_q5_to_float)[0]
        m = len(c)
        A_local = np.zeros((m, m))
        for ii, vi in enumerate(c):
            for jj, vj in enumerate(c):
                if ii == jj:
                    continue
                if qq_distance_sq(verts[vi], verts[vj]) == d1_local:
                    A_local[ii, jj] = 1.0
        L_local.append(np.diag(A_local.sum(axis=1)) - A_local)
    return L_global, L_local


# ----- eigenspaces and zero-extension lift ------------------------------

def _eigenbasis(L: np.ndarray, lam: float, tol: float = 1e-7) -> np.ndarray:
    vals, vecs = np.linalg.eigh(L)
    mask = np.abs(vals - lam) < tol
    return vecs[:, mask]


def _gram_schmidt(M: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    Q: list = []
    for j in range(M.shape[1]):
        v = M[:, j].copy()
        for q in Q:
            v -= (q @ v) * q
        nrm = np.linalg.norm(v)
        if nrm < tol:
            continue
        Q.append(v / nrm)
    return np.column_stack(Q) if Q else np.zeros((M.shape[0], 0))


def _orthogonal_complement(B_full: np.ndarray, B_sub: np.ndarray) -> np.ndarray:
    P = B_full - B_sub @ (B_sub.T @ B_full)
    return _gram_schmidt(P)


def lambda12_lift(state, cosets, L_global, L_local) -> Dict[str, np.ndarray]:
    """Zero-extend local lambda=12 eigenspaces; return E_lifted, E_global(12), E_residual.

    All three bases are returned as orthonormal column matrices.
    """
    n = state["n"]
    cols = []
    for ci, c in enumerate(cosets):
        B = _eigenbasis(L_local[ci], 12.0)
        for col in range(B.shape[1]):
            v = np.zeros(n)
            for k, vidx in enumerate(c):
                v[vidx] = B[k, col]
            cols.append(v)
    B_lifted = _gram_schmidt(np.column_stack(cols))
    B_global = _gram_schmidt(_eigenbasis(L_global, 12.0))
    B_residual = _orthogonal_complement(B_global, B_lifted)
    return {"lifted": B_lifted, "global": B_global, "residual": B_residual}


def lift_residual_at(state, cosets, L_global, L_local, lam: float) -> Dict[str, float]:
    """Compute lift residual + intersection dimensions at an arbitrary lambda.

    Returns:
        lifted_dim:    dim of span of all zero-extended local-lambda eigenvectors.
        global_dim:    dim E_global(lambda).
        intersection_dim: dim ( lift_image  cap  E_global(lambda) ).
        max_residual:  max_i || L_global v_i - lambda v_i ||  for v_i a lifted basis.
                       (Pass-rate metric: equals 0 iff every individual basis vec
                       is a global eigenvector with eigenvalue lambda.)
    """
    n = state["n"]
    cols = []
    for ci, c in enumerate(cosets):
        B = _eigenbasis(L_local[ci], lam)
        for col in range(B.shape[1]):
            v = np.zeros(n)
            for k, vidx in enumerate(c):
                v[vidx] = B[k, col]
            cols.append(v)
    if not cols:
        return {
            "lifted_dim": 0, "global_dim": 0, "intersection_dim": 0,
            "max_residual": float("nan"),
        }
    Bl = _gram_schmidt(np.column_stack(cols))
    Bg = _eigenbasis(L_global, lam)
    if Bg.size > 0:
        Bg = _gram_schmidt(Bg)
        # intersection dim via SVD of B_lifted^T B_global
        M = Bl.T @ Bg
        sigma = np.linalg.svd(M, compute_uv=False) if M.size > 0 else np.array([])
        intersection_dim = int((sigma > 1.0 - 1e-9).sum())
        gdim = int(Bg.shape[1])
    else:
        intersection_dim = 0
        gdim = 0
    res_mat = L_global @ Bl - lam * Bl
    return {
        "lifted_dim": int(Bl.shape[1]),
        "global_dim": gdim,
        "intersection_dim": intersection_dim,
        "max_residual": float(np.linalg.norm(res_mat, axis=0).max()) if Bl.shape[1] > 0 else 0.0,
    }


# ----- exact rational certificate ---------------------------------------

def _integer_adjacency(state, d2_target) -> np.ndarray:
    """Build the exact integer adjacency matrix on V_600 for a fixed d^2 = d2_target."""
    verts = state["verts"]
    n = state["n"]
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if qq_distance_sq(verts[i], verts[j]) == d2_target:
                A[i, j] = 1
    return A


def exact_integer_spectrum_check(state, cosets) -> Dict:
    """Verify, by exact ℚ-rational sympy nullspace, the integer multiplicities
    of the local 24-cell and global 600-cell Laplacians at every integer
    eigenvalue in their respective spectra.

    Returns:
      local_integer_mults: { lambda: int } — coset-wise (same for all 5 cosets).
      global_integer_mults: { lambda: int } — V_600 graph.
    """
    import sympy as sp

    verts = state["verts"]
    n = state["n"]

    # Global integer adjacency
    zero_d2 = (Fraction(0), Fraction(0))
    all_d2: set = set()
    for i in range(n):
        for j in range(n):
            d = qq_distance_sq(verts[i], verts[j])
            if d == zero_d2:
                continue
            all_d2.add(d)
    d1_global = sorted(all_d2, key=_q5_to_float)[0]
    A_global = _integer_adjacency(state, d1_global)
    deg_global = int(A_global.sum(axis=1)[0])
    L_global_int = deg_global * np.eye(n, dtype=np.int64) - A_global

    # Build a local 24-cell Laplacian (same up to isometry for every coset).
    c = cosets[0]
    m = len(c)
    A_local = np.zeros((m, m), dtype=np.int64)
    d_local = (Fraction(1), Fraction(0))
    for ii, vi in enumerate(c):
        for jj, vj in enumerate(c):
            if ii == jj:
                continue
            if qq_distance_sq(verts[vi], verts[vj]) == d_local:
                A_local[ii, jj] = 1
    deg_local = int(A_local.sum(axis=1)[0])
    L_local_int = deg_local * np.eye(m, dtype=np.int64) - A_local

    local_integer_lambdas = (0, 4, 8, 10, 12)
    # Include λ=4,8,10 in the global list too: these are not Laplacian
    # eigenvalues of V_600 but the certificate records the exact nullity 0
    # for use in the selectivity claim.
    global_integer_lambdas = (0, 4, 8, 9, 10, 12, 14, 15)

    local_mults: Dict[int, int] = {}
    Mloc = sp.Matrix(L_local_int.tolist())
    Iloc = sp.eye(m)
    for lam in local_integer_lambdas:
        K = (Mloc - lam * Iloc).nullspace()
        local_mults[lam] = len(K)

    global_mults: Dict[int, int] = {}
    Mglob = sp.Matrix(L_global_int.tolist())
    Iglob = sp.eye(n)
    for lam in global_integer_lambdas:
        K = (Mglob - lam * Iglob).nullspace()
        global_mults[lam] = len(K)

    return {
        "local_integer_mults": local_mults,
        "global_integer_mults": global_mults,
    }


def exact_a5_characters(state, cosets) -> Dict:
    """Exact ℚ-rational A_5 character decomposition for E_lifted, E_residual, E_V600(12).

    Method: every Laplacian is an integer matrix and λ=12 is an integer
    eigenvalue, so all three subspaces have ℚ-rational bases.  For
    every element h of 2I, right-multiplication by h is a permutation
    matrix P_h, and the trace of P_h restricted to a rational subspace
    is itself rational and exactly computable.

    Returns a dict with:
      class_characters: { cls: { 'lifted': int_chi, 'residual': int_chi,
                                  'global': int_chi } } — exact integer characters.
      irrep_multiplicities: { 'lifted': {irrep: int}, ... }
      minus_one_acts_trivially: bool — True iff the trace of the (-1)
        right-multiplication permutation equals the trace of identity
        on each of E_lifted, E_residual, E_V600(12), i.e. iff {±1}
        acts trivially on the λ=12 sector.
    """
    import sympy as sp

    verts = state["verts"]
    idx_of = state["idx_of"]
    n = state["n"]

    # 1) Build exact integer L_global and L_C_i
    zero_d2 = (Fraction(0), Fraction(0))
    all_d2: set = set()
    for i in range(n):
        for j in range(n):
            d = qq_distance_sq(verts[i], verts[j])
            if d == zero_d2:
                continue
            all_d2.add(d)
    sorted_d2 = sorted(all_d2, key=_q5_to_float)
    d1_global = sorted_d2[0]

    A_global = _integer_adjacency(state, d1_global)
    deg_global = int(A_global.sum(axis=1)[0])
    Lg_minus_12 = sp.Matrix(
        (deg_global * np.eye(n, dtype=np.int64) - A_global
         - 12 * np.eye(n, dtype=np.int64)).tolist()
    )

    # 2) Compute exact Q-rational kernel of L_global - 12 I (= E_global(12)).
    #    sympy nullspace gives a list of column vectors with rational entries.
    global_kernel = Lg_minus_12.nullspace()
    # Stack into a 120 x 25 rational matrix B_global
    B_global = sp.Matrix.hstack(*global_kernel) if global_kernel else sp.zeros(n, 0)

    # 3) Compute exact Q-rational E_lifted: for each coset, take rational basis of
    #    E_C(12), zero-extend to R^120.  Stack into 120 x 10 rational matrix.
    lifted_cols = []
    for c in cosets:
        m = len(c)
        A_local = np.zeros((m, m), dtype=np.int64)
        d_local = (Fraction(1), Fraction(0))
        for ii, vi in enumerate(c):
            for jj, vj in enumerate(c):
                if ii == jj:
                    continue
                if qq_distance_sq(verts[vi], verts[vj]) == d_local:
                    A_local[ii, jj] = 1
        deg_local = int(A_local.sum(axis=1)[0])
        Lc_minus_12 = sp.Matrix(
            (deg_local * np.eye(m, dtype=np.int64) - A_local
             - 12 * np.eye(m, dtype=np.int64)).tolist()
        )
        kernel_local = Lc_minus_12.nullspace()
        for vec in kernel_local:
            v_full = sp.zeros(n, 1)
            for k, v_idx in enumerate(c):
                v_full[v_idx, 0] = vec[k, 0]
            lifted_cols.append(v_full)
    B_lifted = sp.Matrix.hstack(*lifted_cols) if lifted_cols else sp.zeros(n, 0)

    # 4) Build exact representation matrices.  Since basis vectors may not be
    #    orthonormal, we use the projection-and-trace formula
    #        chi(h) = trace( (B^T B)^{-1} B^T P_h B )
    #    where P_h is the permutation matrix of right-mult by h.
    def perm_for(h_idx):
        h = verts[h_idx]
        return [idx_of[qq_key(qq_mul(verts[i], h))] for i in range(n)]

    def perm_image_of_basis(B_mat, perm):
        # Apply perm: row i of result = row perm[i] of B_mat.  Equivalently
        # (P_h B)[i,j] = B[perm^{-1}(i), j].  Equivalently row i of P_h B
        # equals row of B at index = inverse permutation at i.
        # Simpler: P_h acts on vectors by (P_h v)[i] = v[perm[i]].
        # So P_h B has row i = B's row perm[i].
        m_rows = B_mat.shape[0]
        new_rows = [B_mat.row(perm[i]) for i in range(m_rows)]
        return sp.Matrix.vstack(*new_rows) if new_rows else B_mat

    # Pre-compute Gram inverse for each basis (rational pseudoinverse via solve)
    def gram_inv(B_mat):
        if B_mat.cols == 0:
            return sp.zeros(0, 0)
        return (B_mat.T * B_mat).inv()

    Gi_lifted = gram_inv(B_lifted)
    Gi_global = gram_inv(B_global)

    # B_residual is the orthogonal complement of B_lifted inside B_global
    # over Q.  Equivalently, columns spanning ker(B_lifted^T) ∩ col(B_global).
    # Use: project B_global to complement of B_lifted (over Q with Gram inverse)
    if B_lifted.cols > 0 and B_global.cols > 0:
        proj = B_lifted * Gi_lifted * (B_lifted.T * B_global)
        B_residual_raw = B_global - proj
        # take a rational column basis via columnspace
        B_residual = sp.Matrix.hstack(*B_residual_raw.columnspace())
    else:
        B_residual = sp.zeros(n, 0)
    Gi_residual = gram_inv(B_residual)

    def chi(h_idx, B_mat, Gi):
        if B_mat.cols == 0:
            return sp.S(0)
        perm = perm_for(h_idx)
        PB = perm_image_of_basis(B_mat, perm)
        M = Gi * (B_mat.T * PB)
        return M.trace()

    # 5) Tally exact characters by A_5 conjugacy class.
    coset_of = {v: ci for ci, c in enumerate(cosets) for v in c}

    def coset_perm(vperm):
        sample = {}
        for v, ci in coset_of.items():
            if ci not in sample:
                sample[ci] = v
        return [coset_of[vperm[sample[ci]]] for ci in range(5)]

    def cyc_struct(perm):
        seen = [False] * len(perm)
        cyc = []
        for i in range(len(perm)):
            if seen[i]:
                continue
            c = 0
            j = i
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                c += 1
            cyc.append(c)
        return tuple(sorted(cyc))

    base_of_struct = {
        (1, 1, 1, 1, 1): "1A", (1, 2, 2): "2A", (1, 1, 3): "3A", (5,): "5*",
    }

    # --- Compute actual A_5 conjugacy classes by orbit enumeration ----
    # 1) Enumerate the 60 distinct 5-letter permutations in the image of 2I -> A_5.
    a5_elts: set = set()
    for h_idx in range(n):
        a5_elts.add(tuple(coset_perm(perm_for(h_idx))))
    assert len(a5_elts) == 60, f"expected 60 A_5 elements, got {len(a5_elts)}"

    def perm_compose(p, q):
        return tuple(p[q[i]] for i in range(len(p)))

    def perm_inverse(p):
        inv = [0] * len(p)
        for i, x in enumerate(p):
            inv[x] = i
        return tuple(inv)

    # 2) Bucket A_5 elements into conjugacy classes via orbit under conjugation.
    a5_class_of_perm: Dict[tuple, str] = {}
    remaining = set(a5_elts)
    class_buckets: Dict[str, set] = {"1A": set(), "2A": set(), "3A": set(),
                                       "5A": set(), "5B": set()}
    # First, assign base classes that are determined uniquely by cycle structure.
    for p in list(remaining):
        struct = cyc_struct(list(p))
        base = base_of_struct[struct]
        if base != "5*":
            class_buckets[base].add(p)
            a5_class_of_perm[p] = base
            remaining.discard(p)
    # Now split the 5-cycles into two A_5-conjugacy classes by orbit enumeration.
    # Pick any unclassified 5-cycle; its A_5-orbit is one class, complement is the other.
    if remaining:
        seed = next(iter(remaining))
        orbit_a = set()
        for g in a5_elts:
            g_inv = perm_inverse(g)
            orbit_a.add(perm_compose(perm_compose(g, seed), g_inv))
        # Sanity: A_5 5-cycle classes have size 12 each.
        if len(orbit_a) != 12:
            raise RuntimeError(
                f"A_5 5-cycle orbit has unexpected size {len(orbit_a)} (expected 12)"
            )
        for p in remaining:
            if p in orbit_a:
                class_buckets["5A"].add(p)
                a5_class_of_perm[p] = "5A"
            else:
                class_buckets["5B"].add(p)
                a5_class_of_perm[p] = "5B"
        # Sanity: class sizes should be 1, 15, 20, 12, 12
        sizes_seen = {cls: len(b) for cls, b in class_buckets.items()}
        expected = {"1A": 1, "2A": 15, "3A": 20, "5A": 12, "5B": 12}
        for cls, want in expected.items():
            if sizes_seen[cls] != want:
                raise RuntimeError(
                    f"A_5 class {cls} has size {sizes_seen[cls]} (expected {want})"
                )

    chi_sum = {cls: {"lifted": sp.S(0), "residual": sp.S(0), "global": sp.S(0)}
               for cls in A5_CLASSES}
    class_count = {cls: 0 for cls in A5_CLASSES}

    minus_one_idx = next(
        i for i in range(n)
        if verts[i] == ((Fraction(-1), Fraction(0)), (Fraction(0), Fraction(0)),
                       (Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))
    )
    plus_one_idx = next(
        i for i in range(n)
        if verts[i] == ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(0)),
                       (Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))
    )

    for h_idx in range(n):
        cperm = tuple(coset_perm(perm_for(h_idx)))
        cls = a5_class_of_perm[cperm]
        cl = chi(h_idx, B_lifted, Gi_lifted)
        cr = chi(h_idx, B_residual, Gi_residual)
        cg = chi(h_idx, B_global, Gi_global)
        chi_sum[cls]["lifted"] += cl
        chi_sum[cls]["residual"] += cr
        chi_sum[cls]["global"] += cg
        class_count[cls] += 1

    class_avg = {cls: {k: (v / class_count[cls] if class_count[cls] > 0 else sp.S(0))
                       for k, v in d.items()}
                 for cls, d in chi_sum.items()}

    # 6) Project against A_5 character table (irrep chars stored as floats in
    #    A5_CHARTAB above; we instead use the exact algebraic form here for 5A/5B).
    #    The two 5-cycle character values are phi = (1+sqrt5)/2 and psi = (1-sqrt5)/2.
    sqrt5 = sp.sqrt(5)
    phi = (1 + sqrt5) / 2
    psi = (1 - sqrt5) / 2
    A5_CHAR_EXACT = {
        "1":  (sp.S(1), sp.S(1), sp.S(1), sp.S(1), sp.S(1)),
        "3":  (sp.S(3), sp.S(-1), sp.S(0), phi, psi),
        "3p": (sp.S(3), sp.S(-1), sp.S(0), psi, phi),
        "4":  (sp.S(4), sp.S(0), sp.S(1), sp.S(-1), sp.S(-1)),
        "5":  (sp.S(5), sp.S(1), sp.S(-1), sp.S(0), sp.S(0)),
    }
    A5_SIZE = sp.S(60)
    decomp = {}
    for name in ("lifted", "residual", "global"):
        decomp[name] = {}
        for irrep, char_row in A5_CHAR_EXACT.items():
            m = sum(sp.S(A5_CLASS_SIZES[i]) * class_avg[A5_CLASSES[i]][name] * char_row[i]
                    for i in range(len(A5_CLASSES))) / A5_SIZE
            m = sp.simplify(m)
            decomp[name][irrep] = m

    # 7) Verify {±1}-triviality on E_V600(12) componentwise:
    #    for every column v of B_global, check (P_{-1} v) - v = 0 exactly in Q^120.
    #    P_{-1} is an involution on R^120, so trace-equality alone implies identity;
    #    we additionally do the strict componentwise check so the proof statement
    #    in the note is supported directly.
    perm_minus = perm_for(minus_one_idx)
    trivial = True
    if B_global.cols > 0:
        PB = perm_image_of_basis(B_global, perm_minus)
        delta = PB - B_global
        for j in range(delta.cols):
            for i in range(delta.rows):
                if delta[i, j] != 0:
                    trivial = False
                    break
            if not trivial:
                break

    # 8) Exact invariance certificate for E_lifted and E_residual.
    #    For each h in 2I, check that P_h B is in the column span of B for both
    #    B = B_lifted and B = B_residual.  Equivalently, check
    #        P_h B = B (B^T B)^{-1} B^T (P_h B)
    #    exactly over Q (no floating point); equality of every coordinate
    #    certifies that P_h preserves the subspace.
    def is_invariant(B_mat, Gi):
        if B_mat.cols == 0:
            return True
        for h_idx in range(n):
            perm_h = perm_for(h_idx)
            PB = perm_image_of_basis(B_mat, perm_h)
            # Projection of PB back into col(B):  B (B^T B)^{-1} B^T PB
            proj = B_mat * (Gi * (B_mat.T * PB))
            delta = PB - proj
            # check componentwise zero over Q
            for j in range(delta.cols):
                for i in range(delta.rows):
                    if delta[i, j] != 0:
                        return False
        return True

    lifted_invariant = is_invariant(B_lifted, Gi_lifted)
    residual_invariant = is_invariant(B_residual, Gi_residual)

    # Cast class averages to int when possible
    def _coerce(x):
        s = sp.simplify(x)
        if s.is_Integer:
            return int(s)
        return float(s)

    class_characters_out = {
        cls: {k: _coerce(v) for k, v in d.items()} for cls, d in class_avg.items()
    }
    irrep_multiplicities_out = {
        name: {irrep: _coerce(m) for irrep, m in decomp[name].items()}
        for name in decomp
    }

    return {
        "class_characters": class_characters_out,
        "irrep_multiplicities": irrep_multiplicities_out,
        "minus_one_acts_trivially": trivial,
        "lifted_2I_invariant": lifted_invariant,
        "residual_2I_invariant": residual_invariant,
    }


def exact_certify_lambda12_lift(state, cosets) -> Dict:
    """Exact rational certificate for the lambda=12 zero-extension lift.

    Both Laplacians have integer entries (the d=1 graphs are 0/1 adjacency
    matrices with integer-valued regular degree).  Hence ker(L - 12 I) is
    a ℚ-rational subspace, and lift exactness is decidable by exact
    rational linear algebra.

    Procedure:
      1. Form the integer matrices L_local_i - 12 I and L_global - 12 I.
      2. Use sympy to compute the ℚ-rational nullspace of each L_local_i - 12 I,
         giving an exact basis for E_{C_i}(12) over ℚ.
      3. Zero-extend each basis vector to a vector in ℚ^120 and apply the
         integer matrix L_global - 12 I.  Each application is exact.
      4. The lift is exactly correct iff the result is the zero vector
         in every coordinate (over ℚ, no rounding).

    Returns a dict with:
      certified: bool — True iff every lifted local-λ=12 eigenvector lies
                 in ker(L_global - 12 I) exactly over ℚ.
      total_local_dim: int — sum_i dim_ℚ E_{C_i}(12) (= 10).
      max_off_residual_components: int — number of nonzero ℚ-components
                 in the (L_global - 12 I) v_lifted vectors, summed over
                 all lifted basis vectors.  Must be 0 for a positive
                 certificate.
    """
    import sympy as sp

    verts = state["verts"]
    n = state["n"]

    # Global d^2 = (3 - sqrt 5)/2 → use the same builder as build_d1_graphs
    zero_d2 = (Fraction(0), Fraction(0))
    all_d2: set = set()
    for i in range(n):
        for j in range(n):
            d = qq_distance_sq(verts[i], verts[j])
            if d == zero_d2:
                continue
            all_d2.add(d)
    sorted_d2 = sorted(all_d2, key=_q5_to_float)
    d1_global = sorted_d2[0]

    A_global = _integer_adjacency(state, d1_global)
    deg_global = int(A_global.sum(axis=1)[0])
    assert all(int(d) == deg_global for d in A_global.sum(axis=1))
    L_global_int = deg_global * np.eye(n, dtype=np.int64) - A_global

    Lg_minus_12 = sp.Matrix((L_global_int - 12 * np.eye(n, dtype=np.int64)).tolist())

    total_local_dim = 0
    nonzero_components = 0
    per_coset = []
    for c in cosets:
        # local d^2 = 1 internally
        m = len(c)
        A_local = np.zeros((m, m), dtype=np.int64)
        d_local = (Fraction(1), Fraction(0))
        for ii, vi in enumerate(c):
            for jj, vj in enumerate(c):
                if ii == jj:
                    continue
                if qq_distance_sq(verts[vi], verts[vj]) == d_local:
                    A_local[ii, jj] = 1
        deg_local = int(A_local.sum(axis=1)[0])
        assert all(int(d) == deg_local for d in A_local.sum(axis=1))
        L_loc_int = deg_local * np.eye(m, dtype=np.int64) - A_local
        Lc_minus_12 = sp.Matrix((L_loc_int - 12 * np.eye(m, dtype=np.int64)).tolist())
        kernel = Lc_minus_12.nullspace()  # list of sympy.Matrix column vectors
        coset_dim = len(kernel)
        total_local_dim += coset_dim
        coset_nonzero = 0
        for vec in kernel:
            # Zero-extend to length-120 sympy column vector
            v_full = sp.zeros(n, 1)
            for k, v_idx in enumerate(c):
                v_full[v_idx, 0] = vec[k, 0]
            applied = Lg_minus_12 * v_full  # exact rational
            for k in range(n):
                if applied[k, 0] != 0:
                    coset_nonzero += 1
        nonzero_components += coset_nonzero
        per_coset.append({"local_kernel_dim": coset_dim, "nonzero_components": coset_nonzero})

    return {
        "certified": (nonzero_components == 0),
        "total_local_dim": total_local_dim,
        "nonzero_components_total": nonzero_components,
        "per_coset": per_coset,
        "method": "sympy ℚ-rational nullspace; integer arithmetic throughout",
    }


# ----- A_5 character decomposition --------------------------------------

PHI = (1.0 + math.sqrt(5.0)) / 2.0
PSI = (1.0 - math.sqrt(5.0)) / 2.0

# A_5 character table (rows = irreps, cols = classes 1A, 2A, 3A, 5A, 5B)
A5_CLASSES = ("1A", "2A", "3A", "5A", "5B")
A5_CLASS_SIZES = (1, 15, 20, 12, 12)
A5_CHARTAB: Dict[str, Tuple[float, ...]] = {
    "1":  (1.0,  1.0,  1.0,  1.0,  1.0),
    "3":  (3.0, -1.0,  0.0,  PHI,  PSI),
    "3p": (3.0, -1.0,  0.0,  PSI,  PHI),
    "4":  (4.0,  0.0,  1.0, -1.0, -1.0),
    "5":  (5.0,  1.0, -1.0,  0.0,  0.0),
}


def _right_mult_perm(state, h_idx: int) -> List[int]:
    verts = state["verts"]
    idx_of = state["idx_of"]
    n = state["n"]
    h = verts[h_idx]
    return [idx_of[qq_key(qq_mul(verts[i], h))] for i in range(n)]


def _cycle_structure(perm: List[int]) -> Tuple[int, ...]:
    n = len(perm)
    seen = [False] * n
    cycles = []
    for i in range(n):
        if seen[i]:
            continue
        c = 0
        j = i
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            c += 1
        cycles.append(c)
    return tuple(sorted(cycles))


def _coset_perm(coset_of: Dict[int, int], vperm: List[int]) -> List[int]:
    """Right-mult on cosets is well-defined; read off using any representative."""
    image = [-1] * 5
    sample = {}
    for v, ci in coset_of.items():
        if ci not in sample:
            sample[ci] = v
    for ci in range(5):
        image[ci] = coset_of[vperm[sample[ci]]]
    return image


def _a5_base_class(struct: Tuple[int, ...]) -> str:
    return {
        (1, 1, 1, 1, 1): "1A",
        (1, 2, 2): "2A",
        (1, 1, 3): "3A",
        (5,): "5*",
    }[struct]


def _character_on(basis: np.ndarray, perm: List[int]) -> float:
    applied = basis[np.asarray(perm), :]
    return float(np.einsum("ij,ij->", basis, applied))


def a5_character_decomposition(
    state,
    cosets,
    bases: Dict[str, np.ndarray],
) -> Dict[str, Dict[str, float]]:
    """Compute A_5 irrep multiplicities for E_lifted, E_residual, E_global(12).

    Method: enumerate all 120 elements of 2I; for each, compute the
    induced permutation of cosets and the trace of its action on each
    basis; tally by A_5 class; project onto A_5 character table.
    """
    n = state["n"]
    coset_of: Dict[int, int] = {v: ci for ci, c in enumerate(cosets) for v in c}

    # Pick reference 5-cycle on coset labels to split 5A / 5B
    five_ref = None
    for h_idx in range(n):
        perm = _right_mult_perm(state, h_idx)
        cperm = _coset_perm(coset_of, perm)
        if _cycle_structure(cperm) == (5,):
            five_ref = cperm
            break
    if five_ref is None:
        raise RuntimeError("no 5-cycle found in coset permutation action")
    t1 = five_ref
    t2 = [t1[i] for i in t1]
    t3 = [t2[i] for i in t1]
    t4 = [t3[i] for i in t1]
    fivea_set = {tuple(t1), tuple(t4)}
    fiveb_set = {tuple(t2), tuple(t3)}

    chi_acc: Dict[str, Dict[str, List[float]]] = {
        cls: {name: [] for name in bases} for cls in A5_CLASSES
    }

    for h_idx in range(n):
        perm = _right_mult_perm(state, h_idx)
        cperm = _coset_perm(coset_of, perm)
        struct = _cycle_structure(cperm)
        base = _a5_base_class(struct)
        if base == "5*":
            tup = tuple(cperm)
            cls = "5A" if tup in fivea_set else ("5B" if tup in fiveb_set else None)
            if cls is None:
                # falls outside both — split arbitrarily but consistently
                cls = "5A"
        else:
            cls = base
        for name, B in bases.items():
            chi_acc[cls][name].append(_character_on(B, perm))

    chi_mean: Dict[str, Dict[str, float]] = {
        cls: {name: (float(np.mean(vs)) if vs else 0.0) for name, vs in d.items()}
        for cls, d in chi_acc.items()
    }

    # Project onto A_5 character table
    A5_SIZE = 60
    decomp: Dict[str, Dict[str, float]] = {name: {} for name in bases}
    for irrep, char_row in A5_CHARTAB.items():
        for name in bases:
            m = sum(
                A5_CLASS_SIZES[i] * chi_mean[A5_CLASSES[i]][name] * char_row[i]
                for i in range(len(A5_CLASSES))
            ) / A5_SIZE
            decomp[name][irrep] = m
    return {"class_characters": chi_mean, "irrep_multiplicities": decomp}


def format_decomposition(mults: Dict[str, float], threshold: float = 0.05) -> str:
    """Render an irrep multiplicity dict as `2*Y5 + ...` if integral, else show raw."""
    label_map = {"1": "Y1", "3": "Y3", "3p": "Y3'", "4": "Y4", "5": "Y5"}
    terms = []
    for irrep_label, m in mults.items():
        if abs(m) < threshold:
            continue
        rounded = int(round(m))
        if abs(m - rounded) < threshold:
            if rounded == 1:
                terms.append(label_map[irrep_label])
            else:
                terms.append(f"{rounded}*{label_map[irrep_label]}")
        else:
            terms.append(f"{m:.3f}*{label_map[irrep_label]}")
    return " + ".join(terms) if terms else "0"


# ----- top-level run -----------------------------------------------------

@dataclass
class KeystoneResult:
    n_vertices: int
    coset_sizes: List[int]
    global_degree: int
    global_edges: int
    local_degrees: List[int]
    local_edges: List[int]
    lambda12_lifted_dim: int
    lambda12_global_dim: int
    lambda12_residual_dim: int
    lambda12_max_residual: float
    negative_controls: List[Dict] = field(default_factory=list)
    irrep_decomposition: Dict[str, Dict[str, float]] = field(default_factory=dict)
    summary_strings: Dict[str, str] = field(default_factory=dict)
    exact_certificate: Dict = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        return (
            "SELECTIVE_SPECTRAL_BRIDGE"
            if self.lambda12_max_residual < 1e-10
            and self.lambda12_lifted_dim == 10
            and self.lambda12_residual_dim == 15
            and self.lambda12_global_dim == 25
            and all(c["max_residual"] > 1.0 for c in self.negative_controls)
            and self.exact_certificate.get("certified", False)
            else "FAIL"
        )

    def to_dict(self) -> Dict:
        return {
            "n_vertices": self.n_vertices,
            "coset_sizes": list(self.coset_sizes),
            "global_degree": self.global_degree,
            "global_edges": self.global_edges,
            "local_degrees": list(self.local_degrees),
            "local_edges": list(self.local_edges),
            "lambda12_lifted_dim": self.lambda12_lifted_dim,
            "lambda12_global_dim": self.lambda12_global_dim,
            "lambda12_residual_dim": self.lambda12_residual_dim,
            "lambda12_max_residual": self.lambda12_max_residual,
            "negative_controls": list(self.negative_controls),
            "irrep_decomposition": self.irrep_decomposition,
            "summary_strings": dict(self.summary_strings),
            "exact_certificate": dict(self.exact_certificate),
            "verdict": self.verdict,
        }


def run_keystone(
    control_lambdas: Tuple[float, ...] = (0.0, 4.0, 8.0, 10.0),
    exact_certificate: bool = True,
) -> KeystoneResult:
    """Run the full keystone protocol and return a structured result.

    If ``exact_certificate`` is True (the default), also runs the ℚ-rational
    sympy-based exactness certificate for the λ=12 lift identity.  Costs
    a few extra seconds but turns the result from "verified to ~1e-15" into
    "verified exactly over ℚ".
    """
    state = build_state()
    cosets = coset_partition(state)
    L_global, L_local = build_d1_graphs(state, cosets)

    # Global graph stats
    deg_global = int(np.diag(L_global)[0])
    edges_global = int(np.diag(L_global).sum() / 2)

    # Local graph stats
    local_degrees = [int(np.diag(L)[0]) for L in L_local]
    local_edges = [int(np.diag(L).sum() / 2) for L in L_local]

    # The lambda=12 lift
    bases = lambda12_lift(state, cosets, L_global, L_local)
    res_max = float(
        np.linalg.norm(L_global @ bases["lifted"] - 12.0 * bases["lifted"], axis=0).max()
    )

    # A_5 character decomposition
    a5 = a5_character_decomposition(state, cosets, bases)
    decomp = a5["irrep_multiplicities"]

    # Negative controls
    controls = []
    for lam in control_lambdas:
        info = lift_residual_at(state, cosets, L_global, L_local, lam)
        info["lambda"] = lam
        controls.append(info)

    summary = {
        "E_lifted": format_decomposition(decomp["lifted"]),
        "E_residual": format_decomposition(decomp["residual"]),
        "E_global": format_decomposition(decomp["global"]),
    }

    if exact_certificate:
        cert: Dict = exact_certify_lambda12_lift(state, cosets)
        cert["integer_spectrum_check"] = exact_integer_spectrum_check(state, cosets)
        cert["a5_exact"] = exact_a5_characters(state, cosets)
    else:
        cert = {}

    return KeystoneResult(
        n_vertices=state["n"],
        coset_sizes=[len(c) for c in cosets],
        global_degree=deg_global,
        global_edges=edges_global,
        local_degrees=local_degrees,
        local_edges=local_edges,
        lambda12_lifted_dim=int(bases["lifted"].shape[1]),
        lambda12_global_dim=int(bases["global"].shape[1]),
        lambda12_residual_dim=int(bases["residual"].shape[1]),
        lambda12_max_residual=res_max,
        negative_controls=controls,
        irrep_decomposition=decomp,
        summary_strings=summary,
        exact_certificate=cert,
    )
