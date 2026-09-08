orders := [];

for Flist in Fs do
    F := Matrix(Integers(), Flist);
    L := LatticeWithGram(F);
    G := AutomorphismGroup(L);
    Append(~orders, #G);
end for;

orders;

