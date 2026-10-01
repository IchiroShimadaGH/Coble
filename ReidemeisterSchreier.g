#Read("ReidemeisterSchreier.g");

# Made by Codex


# Reidemeister-Schreier generators for the stabilizer of 1.
# Input: [p1,...,pN], images of the ORIGINAL generators [g1,...,gN].
# GAP right-action convention: x^(p*q) = (x^p)^q.
# No presentation of G and no faithful action are required.

ReidemeisterSchreierStabilizer := function(perms)
  #
  local N, F, letters, orbit, position, transversal, head, x, i, y, j, 
  words, w, integerWords, rrec, tyinv, tj;
  #
  if not IsList(perms) or not ForAll(perms, IsPerm) then
    Error("Expected a list [p1,...,pN] of permutations");
  fi;
  #
  N:=Length(perms); 
  F:=FreeGroup(N); 
  letters:=GeneratorsOfGroup(F);
  #
  orbit:=[1]; 
  position:=[]; 
  position[1]:=1;
  transversal:=[One(F)]; 
  #
  head:=1;
  # Breadth-first Schreier tree. Positive generators suffice because the
  # permutation orbit is finite. Each t_x satisfies 1^image(t_x)=x.
  while head<=Length(orbit) do
    x:=orbit[head];
    for i in [1..N] do
    y:=x^perms[i];
    if not IsBound(position[y]) then
      Add(orbit, y); position[y]:=Length(orbit);
      Add(transversal, transversal[head]*letters[i]);
    fi;
    od;
    head:=head+1;
  od;
  #
  words:=[];
  for j in [1..Length(orbit)] do
    x:=orbit[j];
    tj:=transversal[j];
    for i in [1..N] do
    y:=x^perms[i];
    tyinv:=transversal[position[y]]^(-1);
    w:=tj*letters[i]*tyinv;
    # Omit only freely trivial words. In particular,  NEVER omit a word
    # because its permutation image is identity: it may be a kernel element.
    if w<>One(F) then Add(words, w); fi;
    od;
  od;
  integerWords:=List(words, LetterRepAssocWord);
  #
  rrec:=rec(
      freeGroup:=F, 
      generators:=letters, 
      words:=words,
      integerWords:=integerWords, 
      orbit:=orbit,
      transversal:=transversal, 
      index:=Length(orbit));
      #
  return (rrec);
end;

# Minimal interface: returns just the list of words in f1,...,fN.
RSStabilizerWords := function(perms)
 return ReidemeisterSchreierStabilizer(perms).words;
end;
