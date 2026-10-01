
#Read("Task20260930Ub.g");

Read("NonSingOverLats.g");

readdata("FoundShRecs0930Ua");

Printn("FoundShRecs0930Ua", Length(FoundShRecs0930Ua));

FoundShRecsOverLats:=[];

counter:=0;

for trec in FoundShRecs0930Ua do
  counter:=counter+1; 
  Printn("_____________");
  OLrecs:=NonSingOverLats(trec.GramS, trec.h, trec.AutShgens);
  trec.OLrecs:=OLrecs;
  Printn(counter, Length(OLrecs), List(OLrecs, xx->xx.extdeg));
od;
##########

