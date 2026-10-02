Read("outputs/random_rank19.g");;
a := R19Search(rec(target:=100,maxTrials:=1000,file:="work/chain.g",progress:=0));;
if a.new <> 100 then Error("chain count"); fi;
for e in R19Saved do
 if not R19Cyclic(e.gram) then Error("chain cyclic"); fi;
 if DeterminantIntMat(e.gram) <> e.determinant then Error("det"); fi;
od;
b := R19Search(rec(mode:="dense",target:=30,maxTrials:=1000,file:="work/dense.g",progress:=0));;
if b.new <> 30 then Error("dense count"); fi;
for e in R19Saved do
 for i in [1..19] do
  if e.gram[i][i] mod 2 <> 0 or e.gram[i][i]-Sum([1..19],j->AbsInt(e.gram[i][j]))+e.gram[i][i] < 4 then Error("bound"); fi;
 od;
 if not R19Cyclic(e.gram) then Error("dense cyclic"); fi;
od;
c := R19Search(rec(target:=10,file:="work/chain.g",resume:=true,seed:=123,progress:=0));;
if c.total <> 110 then Error("resume"); fi;
Print("TESTS PASSED\n");
QUIT;
