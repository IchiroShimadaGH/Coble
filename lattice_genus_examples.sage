# Run in this directory: sage lattice_genus_examples.sage
# Also usable from a Sage notebook with this directory on Python's sys.path.
from lattice_genus_backend import even_genus_representatives

grams = even_genus_representatives([2,17], [30], matrix(QQ, [[17/30]]))
assert len(grams) == 1
print("signature (2,17):", len(grams), "class")
assert grams[0].nrows() == 19 and grams[0].det() == -30

grams = even_genus_representatives([1,1], [229], [[-2/229]])
assert len(grams) == 2
print("binary example:", len(grams), "classes")
for gram in grams:
    print(gram)

grams = even_genus_representatives([1,2], [30], [[17/30]])
assert grams == []
print("empty genus:", len(grams), "classes")

p = 1513
q = diagonal_matrix(QQ, [1/2, 1/(2*p), -1/(2*p^2)])
grams = even_genus_representatives([2,1], [2,2*p,2*p^2], q)
assert len(grams) == 4
print("ternary example:", len(grams), "classes")
