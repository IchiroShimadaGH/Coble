Read("outputs/random_rank19.g");;
r := R19Search(rec(target:=5,diagonalMax:=6,file:="work/duplicates_results.g",progress:=0));;
if r.total <> 5 or Length(Set(List(R19Saved,e->e.gram))) <> 1 then Error("duplicates test"); fi;
Print("DUPLICATES TEST PASSED\n");
QUIT;
