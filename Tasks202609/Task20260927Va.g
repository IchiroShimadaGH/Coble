
#Read("Task20260927Va.g");


Read("MinimalWgs.g");



MinimalAbsWgsK3Recs:=[[rec(wg:=[], stab:=[()])]];

for newkk in [2..15] do 
  oldminabswgs:=
  List(MinimalAbsWgsK3Recs[Length(MinimalAbsWgsK3Recs)], xx->xx.wg);
  newminabswgrecs:= Enhanced2MakeMinimalAbsWgs(oldminabswgs, newkk);
  Printn(newkk, Length(newminabswgrecs));
  Add(MinimalAbsWgsK3Recs, newminabswgrecs);
  savedata(MinimalAbsWgsK3Recs);
od;

##########

