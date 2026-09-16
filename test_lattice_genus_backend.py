"""Run with: sage -python test_lattice_genus_backend.py"""
from sage.all import ZZ, QQ, matrix, diagonal_matrix, BinaryQF_reduced_representatives
from sage.quadratic_forms.genera.genus import Genus
from lattice_genus_backend import even_lattice_genus, binary_representatives


def check(signature, orders, form, count, name):
    result = even_lattice_genus(signature, orders, form)
    assert result["count"] == count, (name, result)
    print(name, count, result["method"], result["spinorPrimes"], flush=True)
    return result


check([2, 17], [30], [[QQ(17)/30]], 1, "previous rank-19 example")
check([2, 15], [3, 3, 42],
      [[-QQ(4)/3, 0, -QQ(2)/3], [0, 0, -QQ(1)/3],
       [-QQ(2)/3, -QQ(1)/3, -QQ(25)/42]], 1, "previous rank-17 example")
check([1, 2], [30], [[QQ(17)/30]], 0, "incompatible signature")
check([2, 1], [], [], 0, "impossible unimodular signature")
check([1, 1], [], [], 1, "hyperbolic plane")
check([1, 1], [2, 2], [[0, QQ(1)/2], [QQ(1)/2, 0]], 1, "U(2), square discriminant")
check([1, 1], [229], [[-QQ(2)/229]], 2, "binary genus with two GL classes")
for p, count in [(17, 2), (1513, 4)]:
    result = check([2, 1], [2, 2*p, 2*p*p],
                  [[QQ(1)/2, 0, 0], [0, QQ(1)/(2*p), 0],
                   [0, 0, -QQ(1)/(2*p*p)]], count, "ternary spinor test p="+str(p))
    expected_genus = Genus(diagonal_matrix(ZZ, [2, 2*p, -2*p*p]))
    assert all(Genus(matrix(ZZ, a)) == expected_genus for a in result["grams"])

# Compare our integer-bound binary enumeration against Sage's independent
# public enumeration for a range of small discriminants, including squares.
for d in range(1, 81):
    if d % 4 not in (0, 1):
        continue
    public = BinaryQF_reduced_representatives(d, primitive_only=False, proper=False)
    genera = []
    for f in public:
        G = Genus(matrix(ZZ, [[2*f[0], f[1]], [f[1], 2*f[2]]]))
        if G not in genera:
            genera.append(G)
    for G in genera:
        expected = sum(Genus(matrix(ZZ, [[2*f[0], f[1]], [f[1], 2*f[2]]])) == G
                       for f in public)
        assert len(binary_representatives(G)) == expected, d
print("binary cross-checks through discriminant 80 passed", flush=True)

for sig, inv, form in [([2, 1], [3], [[QQ(1)/3]]),
                       ([2, 1], [6], [[QQ(2)/3]]),
                       ([2, 1], [3], [[QQ(1)/2]])]:
    try:
        even_lattice_genus(sig, inv, form)
    except ValueError:
        pass
    else:
        raise AssertionError("invalid discriminant input accepted")
print("invalid inputs rejected; all tests passed", flush=True)
