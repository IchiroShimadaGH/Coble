#Read("Task20260904Ub.g");

Read("Make_nSaSvSss.g");
Read("Adj.g");

readdata("Borrec");
readdata("nSaSvSss0");
readdata("adjrecs0");



numdata0:=Collected(List(nSaSvSss0, trec->[trec.nS, trec.aS,  Length(trec.vSliftss)]));
Printn(numdata0);
Printn("___________");
nSaSvSrecs1:=[];

numdatas:=[];

for tadjrec in adjrecs0 do 
  tweyl:=tadjrec.newweyl;
  newsSaSvSss:=Make_nSaSvSss(Borrec, tweyl);
  check_nSaSvSss(Borrec, tweyl, newsSaSvSss);
  numdata:=Collected(List(newsSaSvSss, trec->[trec.nS, trec.aS,  Length(trec.vSliftss)]));
  Printn(numdata);
  Add(numdatas,numdata);
  Add(nSaSvSrecs1,  newsSaSvSss);
od;

Printn(Collected(numdatas));
savedatap(nSaSvSrecs1);


####################