
#Read("Task20260927Ub.g");


Read("MinimalWgs.g");



MinimalAbsWgsRecs:=[[rec(wg:=[], stab:=[()])]];

for newkk in [2..15] do 
  oldminabswgs:=List(MinimalAbsWgsRecs[Length(MinimalAbsWgsRecs)], xx->xx.wg);
  newminabswgrecs:= 
EnhancedMakeMinimalAbsWgs(oldminabswgs, newkk);
  Printn(newkk, Length(newminabswgrecs));
  Add(MinimalAbsWgsRecs, newminabswgrecs);
  savedata(MinimalAbsWgsRecs);
od;

##########

