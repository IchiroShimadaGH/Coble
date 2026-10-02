# GAP 4; no packages needed. Norm means x * G * x^t.
# Read this file, then call R19Search(options). See README_ja.md.

if not IsBound(R19Saved) then R19Saved := []; fi;

R19Candidate := function(rng, mode, diagonalMax, offBound)
    local n, g, i, j, a, s, lo;
    n := 19;
    g := List([1..n], i -> List([1..n], j -> 0));
    if mode = "chain" then
        for i in [1..n-1] do
            g[i][i+1] := 1; g[i+1][i] := 1;
        od;
        for i in [1..n] do
            g[i][i] := 2 * Random(rng, [3..QuoInt(diagonalMax,2)]);
        od;
    elif mode = "dense" then
        for i in [1..n-1] do
            for j in [i+1..n] do
                a := Random(rng, [-offBound..offBound]);
                g[i][j] := a; g[j][i] := a;
            od;
        od;
        for i in [1..n] do
            s := Sum([1..n], j -> AbsInt(g[i][j]));
            lo := QuoInt(s+5,2); # ceil((s+4)/2)
            # diagonalMax bounds the extra even diagonal slack in dense mode.
            g[i][i] := 2*(lo + Random(rng,[0..QuoInt(diagonalMax-6,2)]));
        od;
    else
        Error("mode must be chain or dense");
    fi;
    return g;
end;

R19Cyclic := function(g)
    local s;
    s := SmithNormalFormIntegerMat(g);
    return ForAll([1..18], i -> AbsInt(s[i][i]) = 1);
end;

R19Search := function(opt)
    local cfg, name, rng, seen, tries, accepted, g, det, entry, start;
    cfg := rec(mode := "chain", target := 1000, maxTrials := 100000,
               distinctDet := false, seed := 20261002, diagonalMax := 40, offBound := 1,
               file := "rank19_results.g", resume := false, progress := 1000);
    for name in RecNames(opt) do
        if not IsBound(cfg.(name)) then Error("Unknown option: ",name); fi;
        cfg.(name) := opt.(name);
    od;
    if not cfg.mode in ["chain","dense"] or
       not IsInt(cfg.target) or cfg.target < 1 or
       not IsInt(cfg.maxTrials) or cfg.maxTrials < 1 or
       not IsInt(cfg.diagonalMax) or cfg.diagonalMax < 6 or
       not IsInt(cfg.offBound) or cfg.offBound < 0 or
       not IsInt(cfg.progress) or cfg.progress < 0 then
        Error("Invalid search options");
    fi;
    if IsExistingFile(cfg.file) then
        if not cfg.resume then
            Error("Output exists; choose a new file or use resume := true");
        fi;
        R19Saved := [];
        Read(cfg.file);
    else
        R19Saved := [];
        PrintTo(cfg.file, "# Rank 19 even positive definite root-free lattices\n",
                "R19Saved := [];;\n");
    fi;
    seen := Set(List(R19Saved, e -> e.determinant));
    rng := RandomSource(IsMersenneTwister, cfg.seed);
    tries := 0; accepted := 0; start := Runtime();
    while accepted < cfg.target and tries < cfg.maxTrials do
        tries := tries+1;
        g := R19Candidate(rng,cfg.mode,cfg.diagonalMax,cfg.offBound);
        det := DeterminantIntMat(g);
        if (not cfg.distinctDet or not det in seen) and (cfg.mode = "chain" or R19Cyclic(g)) then
            entry := rec(gram := g, determinant := det,
                         discriminantInvariantFactors :=
                             Concatenation(List([1..18],i->1),[det]),
                         normLowerBound := 4, mode := cfg.mode,
                         seed := cfg.seed, trial := tries);
            AppendTo(cfg.file, "Add(R19Saved, ", entry, ");;\n");
            Add(R19Saved,entry); AddSet(seen,det);
            accepted := accepted+1;
        fi;
        if cfg.progress > 0 and tries mod cfg.progress = 0 then
            Print("trials=",tries," new=",accepted,
                  " total=",Length(R19Saved)," elapsed_ms=",Runtime()-start,"\n");
        fi;
    od;
    Print("Finished: trials=",tries," new=",accepted,
          " total=",Length(R19Saved)," elapsed_ms=",Runtime()-start,"\n");
    return rec(trials := tries, new := accepted, total := Length(R19Saved),
               targetReached := accepted = cfg.target, file := cfg.file);
end;
