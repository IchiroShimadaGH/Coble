#Read("Task20260909Ue.g");


orders:=
[ 8561413324800, 1712282664960, 743178240, 1698693120, 127401984, 99532800,
298598400, 5096079360, 743178240, 31708938240, 452984832, 28665446400,
16307453952, 3483648000, 2090188800, 707788800, 67947724800, 67947724800,
133772083200, 9555148800, 535088332800, 570760888320, 321052999680,
2853804441600, 64210599936000 ]
; #made by load "Task20260909Uf.m"; 

Read("MassFormula.txt");

if Set(Fs1, Mass)<>[Sum(orders, xx->1/xx)] then beep(41764761); fi;

Rrecs:=[];
for pos in [1..25] do 
  trec:=easybasisrecs1[pos];
  if Fs1[pos]<>trec.Gram then beep(881121); fi;
  tGram:=trec.newGram;
  ogsize:=orders[pos];
  keretasize:=KerOLtoOqL(tGram).size;
  imetaindex:=ogsize/keretasize;
  if not IsInt(imetaindex) then beep(716251); fi;
  ADEtype:=trec.ADEtype;
  svsrec:=ShortestVectors(tGram, 4);
  m2:=2*OccurNumb(svsrec.norms, 2);
  m4:=2*OccurNumb(svsrec.norms, 4);
  Rrec:=rec(
    Gram:=tGram, 
    ogsize:=ogsize,
    imetaindex:=imetaindex, 
    ADEtype:= ADEtype, 
    ms:=[m2, m4], 
  );
  Printn(Collected(ADEtype), ogsize, imetaindex);
  Add(Rrecs, Rrec);
od;

savedata(Rrecs);

Printn(Minimum(List(Rrecs, xx->xx.imetaindex)));