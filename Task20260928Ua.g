
#Read("Task20260928Ua.g");


Read("MinimalWgs.g");



MinAbsWgsK3Recs:=[[rec(wg:=[], stab:=[()])]];

for newkk in [2..15] do 
  oldminabswgs:=
  List(MinAbsWgsK3Recs[Length(MinAbsWgsK3Recs)], xx->xx.wg);
  newminabswgrecs:= Enhanced2MakeMinimalAbsWgs(oldminabswgs, newkk);
  Printn(newkk, Length(newminabswgrecs));
  Add(MinAbsWgsK3Recs, newminabswgrecs);
  savedata(MinAbsWgsK3Recs);
od;

##########

