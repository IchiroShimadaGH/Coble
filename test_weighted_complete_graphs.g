Read("weighted_complete_graphs.g");

# Independent orbit count using Burnside's lemma.
# A cycle of length d in the edge permutation has 5 fixed choices ordinarily.
# After global complementation it has 5 choices if d is even, and 1 if d is odd.
WCGCountByBurnside := function(k)
    local edges, total, p, ep, cycles;
    edges := Combinations([1..k],2);
    total := 0;
    for p in SymmetricGroup(k) do
        ep := PermList(List(edges,e -> Position(edges,Set(e,i -> i^p))));
        cycles := Cycles(ep,[1..Length(edges)]);
        total := total + 5^Length(cycles)
                       + 5^Number(cycles,c -> Length(c) mod 2=0);
    od;
    return total/(2*Factorial(k));
end;

WCGRunTests := function()
    local k, reps, edges, maps, brute, w, map, image, best, n;
    for k in [0..5] do
        reps := WeightedCompleteGraphClasses(k);
        n := WCGCountByBurnside(k);
        if Length(reps)<>n or Length(Set(reps))<>n then
            Error("Incorrect representative count at k=",k);
        fi;
        # Exhaustive comparison to the full domain, independently of prefix pruning.
        if k<=4 then
            edges := Combinations([1..k],2);
            maps := List(Elements(SymmetricGroup(k)),p ->
                List(edges,e -> Position(edges,Set(e,i -> i^p))));
            brute := [];
            for w in Tuples([0..4],Length(edges)) do
                best := w;
                for map in maps do
                    image := w{map};
                    best := Minimum(best,image,List(image,a -> 4-a));
                od;
                AddSet(brute,best);
            od;
            if reps<>brute then Error("Exhaustive comparison failed"); fi;
        fi;
        Print("k=",k,": ",n," classes; passed\n");
    od;
    if WeightedCompleteGraphMatrix(3,[0,1,4])<>[[0,0,1],[0,0,4],[1,4,0]] then
        Error("Matrix conversion failed");
    fi;
    Print("All tests passed. k=6 count (no enumeration): ",WCGCountByBurnside(6),"\n");
end;
WCGRunTests();
