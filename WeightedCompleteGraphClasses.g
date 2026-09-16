#Read("WeightedCompleteGraphClasses.g");

# GAP 4, no additional packages required.
# Representatives are lexicographically minimal edge-weight lists.
# Edge order: Combinations([1..k],2).
# Equivalence: vertex permutations, together with ONE global w -> 4-w.

WeightedCompleteGraphClasses := function(k)
 local edges, m, maps, p, w, reps, CanStillBeCanonical, Search;
 if not IsInt(k) or k < 0 then
  Error("k must be a nonnegative integer");
 fi;
 edges := Combinations([1..k],2);
 m := Length(edges);
 maps := [];
 for p in SymmetricGroup(k) do
  AddSet(maps,List(edges,e -> Position(edges,Set(e,i -> i^p))));
 od;
 w := [];
 reps := [];
 #
 CanStillBeCanonical := function(d)
  local map, flip, i, j, a;
  for map in maps do
for flip in [false,true] do
 for i in [1..d] do
  j := map[i];
  # An unknown earlier entry prevents lexicographic comparison.
  if j > d then break; fi;
  a := w[j];
  if flip then a := 4-a; fi;
  if w[i] < a then break; fi;
  if w[i] > a then return false; fi;
 od;
od;
  od;
  return true;
 end;

 Search := function(d)
  local a;
  if d = m then
Add(reps,ShallowCopy(w));
return;
  fi;
  for a in [0..4] do
w[d+1] := a;
if CanStillBeCanonical(d+1) then Search(d+1); fi;
  od;
  Unbind(w[d+1]);
 end;
 Search(0);
 return reps;
end;

# Convert a representative into Gram matrix.

WeightedCompleteGraphToGram := function(k,w)
  #
 local edges, Gram, t, i, j, ii;
 if not IsInt(k) or k < 0 then
  Error("k must be a nonnegative integer");
 fi;
 edges := Combinations([1..k],2);
 if not IsList(w) or Length(w) <> Length(edges)
 or not ForAll(w,x -> IsInt(x) and x >= 0 and x <= 4) then
  Error("w must contain one weight in [0..4] for each edge");
 fi;
 Gram :=(-2)*IdentityMat(k+1);
 Gram[1][1]:=2;
 for ii in [2..k+1] do 
  Gram[1][ii] := 2; 
  Gram[ii][1] := 2;
 od;
 for t in [1..Length(edges)] do
  i := edges[t][1]; j := edges[t][2];
  Gram[i+1][j+1] := w[t]; 
  Gram[j+1][i+1] := w[t];
 od;
 return Gram;
end;
