"""SageMath backend for even_lattice_genus.g (SageMath 10.9 tested).

Quadratic discriminant values are in Q/2Z; matrices are bilinear Gram
matrices, NOT half-Gram matrices. Ordinary integral isometry is used.
"""
import json
import sys
from pathlib import Path

from sage.all import (ZZ, QQ, Integer, matrix, prod, QuadraticForm, BinaryQF,
                      next_prime)
from sage.modules.torsion_quadratic_module import TorsionQuadraticForm
from sage.quadratic_forms.genera.genus import Genus


def discriminant_input(signature, orders, entries):
    if len(signature) != 2 or any(not isinstance(x, (int, Integer)) or x < 1 for x in signature):
        raise ValueError("signature must be [positive_rank, negative_rank], both positive")
    if any(not isinstance(d, (int, Integer)) or d < 2 for d in orders):
        raise ValueError("cyclic orders must be integers >= 2; use [] for the trivial group")
    m = len(orders)
    if hasattr(entries, "rows"):
        entries = entries.rows()
    if len(entries) != m or any(len(row) != m for row in entries):
        raise ValueError("form dimensions must match the number of cyclic factors")
    q = matrix(QQ, m, m, [QQ(x) for row in entries for x in row])
    if q != q.transpose():
        raise ValueError("the discriminant matrix must be symmetric")
    for i, d in enumerate(orders):
        if any(d*q[i, j] not in ZZ for j in range(m)):
            raise ValueError("bilinear form is incompatible with the stated cyclic orders")
        if d*d*q[i, i] not in 2*ZZ:
            raise ValueError("quadratic form is not well defined modulo 2Z")
    D = TorsionQuadraticForm(q)
    if D.cardinality() != prod(orders):
        raise ValueError("the bilinear discriminant form is degenerate on the stated group")
    # is_genus also checks the quadratic value modulus.
    return D


def binary_representatives(G):
    """Exact reduced enumeration, including square discriminants and content>1.

    Work with ax^2+bxy+cy^2, whose lattice Gram matrix is [[2a,b],[b,2c]].
    Integer square roots avoid floating-point bounds in reduced enumeration.
    """
    d = ZZ(-G.determinant())
    if d % 4 not in (0, 1):
        raise ArithmeticError("an even binary lattice must have discriminant 0 or 1 mod 4")
    root = d.isqrt()
    if root*root == d:
        candidates = ((a, root, ZZ(0)) for a in range(int((-root)//2+1), int(root//2+1)))
    else:
        def generate():
            for b0 in range(1, int(root)+1):
                b = ZZ(b0)
                if (d-b*b) % 4:
                    continue
                ac = (d-b*b)//4
                low = (root-b)//2+1
                high = ac.isqrt()
                for a in ac.divisors():
                    if not low <= a <= high:
                        continue
                    c = -ac//a
                    yield a, b, c
                    yield -a, b, -c
                    yield c, b, a
                    yield -c, b, -a
        candidates = generate()
    forms, grams = [], []
    for a, b, c in candidates:
        M = matrix(ZZ, [[2*a, b], [b, 2*c]])
        if Genus(M) != G:
            continue
        f = BinaryQF([a, b, c])
        if any(f.is_equivalent(g, proper=False) for g in forms):
            continue
        forms.append(f)
        grams.append(M)
    if not grams:
        raise ArithmeticError("binary enumeration produced no representative for a nonempty genus")
    return grams


def spinor_representatives(G):
    """One p-neighbor for EACH nonidentity improper spinor class.

    For indefinite rank >=3, improper spinor classes equal ordinary isometry
    classes. We enumerate every quotient element, not merely generators.
    Sage 10.9's spinor_generators/representatives must not be used for this:
    its subgroup update includes all quotient generators and may stop early.
    """
    A, K = G._improper_spinor_kernel()
    quotient = A.quotient(K)
    count = ZZ(quotient.order())
    seen = [quotient.one()]
    grams = [G.representative()]
    primes = [0]
    base = QuadraticForm(ZZ, grams[0])
    p = ZZ(2)
    while len(seen) < count:
        p = next_prime(p)
        if G.determinant() % p == 0:
            continue
        coset = quotient(A.delta(p))
        if coset in seen:
            continue
        v = base.find_primitive_p_divisible_vector__next(p)
        if v is None:
            raise ArithmeticError("no isotropic vector at a good odd prime")
        neighbor = base.find_p_neighbor_from_vec(p, v)
        grams.append(neighbor.Hessian_matrix().change_ring(ZZ))
        seen.append(coset)
        primes.append(int(p))
    if len(grams) != count:
        raise ArithmeticError("spinor quotient and representative counts disagree")
    return grams, primes


def even_lattice_genus(signature, orders, entries):
    D = discriminant_input(signature, orders, entries)
    signature = tuple(signature)
    if not D.is_genus(signature):
        return dict(count=0, grams=[], method="existence criterion", spinorPrimes=[])
    G = D.genus(signature)
    if not G.is_even() or G.signature_pair() != signature:
        raise ArithmeticError("constructed genus has incorrect parity or signature")
    target = D.normal_form().gram_matrix_quadratic()
    if G.discriminant_form().normal_form().gram_matrix_quadratic() != target:
        raise ArithmeticError("constructed genus has incorrect discriminant form")
    n = sum(signature)
    primes = []
    if n >= len(D.invariants()) + 2:
        grams = [G.representative()]
        method = "Nikulin uniqueness"
    elif n == 2:
        grams = binary_representatives(G)
        method = "binary reduced forms, GL(2,Z) equivalence"
    else:
        grams, primes = spinor_representatives(G)
        method = "all improper spinor classes"
    for M in grams:
        if M.nrows() != n or M != M.transpose() or any(a % 2 for a in M.diagonal()):
            raise ArithmeticError("invalid even Gram matrix")
        if M.det() != (-1)**signature[1]*prod(orders) or Genus(M) != G:
            raise ArithmeticError("returned matrix is not in the requested genus")
    return dict(count=len(grams), grams=[[[int(x) for x in row] for row in M.rows()]
                                       for M in grams], method=method, spinorPrimes=primes)


def even_genus_representatives(signature, orders, form):
    """Sage-native interface: return all representative ZZ Gram matrices.

    Return [] precisely when no lattice exists. Accepts Sage integers and
    either a Sage rational matrix or a nested list. Invalid data raise errors.
    The class number is len(the returned list).
    """
    answer = even_lattice_genus(signature, orders, form)
    return [matrix(ZZ, gram) for gram in answer["grams"]]


def gap_result(result):
    return ("return rec(count := %d, grams := %s, method := %s, spinorPrimes := %s);\n"
            % (result["count"], json.dumps(result["grams"]),
               json.dumps(result["method"]), json.dumps(result["spinorPrimes"])))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: sage -python lattice_genus_backend.py input.json output.g")
    data = json.loads(Path(sys.argv[1]).read_text())
    answer = even_lattice_genus(data["signature"], data["invariants"], data["form"])
    Path(sys.argv[2]).write_text(gap_result(answer))
