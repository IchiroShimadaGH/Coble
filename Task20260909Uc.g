#Read("Task20260909Uc.g");


#
#check NewOGLat and FindIsomBasisRecs
#

Read("NewOGLat.g");

rtypes:=["A1", "A1","A1","A1","A2", "A3","A1", "A2", "A3",
"A1", "A2", "A3", "A4","A5", "A6", "A7", "A8", "A9", "A10", 
"D4", "D5", "D6", "D7", "D8", "D9", "D10", "E6", "E7", "E8"];


for ttt in [1..10000] do
  ttype:=[];
  trank:=0;
  while trank<10 do 
    Add(ttype, Random(rtypes));
    trank:=ADEHrank(ttype);
  od;
  if trank>16 then continue; fi;
  ttGram:=ADEHGram(ttype);
  tU:=RandomUnimodMat(trank);
  tGram:=TMTTmult(tU, ttGram);
  Printn(Collected(ttype));
  st:=Runtime();
  tbasisrec:=BasisRec(tGram, 30);
  tOGrec:=NewOGLat(tbasisrec);
  tsize:=tOGrec.size;
  if tsize<>WeylGroupOrder(ttype)*ADEDiagramAutOrder(ttype) then beep(6276278); fi;
  Printn("_____________", Collected(ttype), tsize, TimeToString(Runtime()-st));
od;




#####





