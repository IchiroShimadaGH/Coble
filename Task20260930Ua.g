
#Read("Task20260930Ua.g");


Read("MinimalWgs.g");



MinAbsWgsK3Recs:=[[rec(wg:=[], stab:=[()])]];
FoundShRecs:=[];

for newkk in [2..8] do 
  oldminabswgs:=
  List(MinAbsWgsK3Recs[Length(MinAbsWgsK3Recs)], xx->xx.wg);
  newminabswgrecs:=EnhancedMinimalWgs(oldminabswgs, newkk);
  Printn(newkk, Length(newminabswgrecs));
  Add(MinAbsWgsK3Recs, newminabswgrecs);
  savedataas(FoundShRecs, "FoundShRecs0930Ua");
  savedataas(MinAbsWgsK3Recs, "MinAbsWgsK3Recs0930Ua");
od;

##########

