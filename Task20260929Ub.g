
#Read("Task20260929Ub.g");


Read("MinimalWgs.g");



MinAbsWgsK3Recs:=[[rec(wg:=[], stab:=[()])]];
FoundShRecs:=[];

for newkk in [2..8] do 
  oldminabswgs:=
  List(MinAbsWgsK3Recs[Length(MinAbsWgsK3Recs)], xx->xx.wg);
  newminabswgrecs:=NewEnhancedMinimalWgs(oldminabswgs, newkk);
  Printn(newkk, Length(newminabswgrecs));
  Add(MinAbsWgsK3Recs, newminabswgrecs);
  savedata(FoundShRecs);
  savedata(MinAbsWgsK3Recs);
od;

##########

