#Read("Task20260918Ub.g");


nn := 9;
hh := Concatenation([3], List([1..8], i -> -1));


codedim:=0;
for codes in CodesByDimension do 
  ttypes:=[];
  for code in codes do 
    latrec:=OverlatticeGramMatrix(code);
    thh:=hh*InverseMat(latrec.basis);
    if not ForAll(thh, IsInt) then buzz(51231); fi;
    tGram:=latrec.Gram;
    if thh*tGram*thh <>2 then buzz(816521); fi;
    singrs:=AffESstd(tGram, thh, 0, -2, true);
    if singrs=[] then 
      Printn(code);
    else 
      Add(ttypes, GetRootType(tGram, singrs));
    fi;
  od;
  Printn("codedim", codedim, "____ done", Collected(ttypes));
  codedim:=codedim+1;
od;

# gap> Read("Task20260918Ub.g");
# [ 0 ] 
# codedim 0 ____ done [  ] 
# codedim 1 ____ done [ [ [ "A1", "A1", "A1", "A1", "A1", "A1", "A1", "A1" ], 135 ] ] 
# codedim 2 ____ done [ [ [ "D4", "D4" ], 1575 ] ] 
# codedim 3 ____ done [ [ [ "D8" ], 2025 ] ] 
# codedim 4 ____ done [ [ [ "E8" ], 270 ] ] 
