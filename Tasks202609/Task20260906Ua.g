#Read("Task20260906Ua.g");

#
#check NewOGLat and FindIsomBasisRecs
#

Read("NewOGLat.g");

for ttt in [1..10000] do
  for tdim in [2..5] do
    tGram:=RandomPosLat(tdim, [-2,-1,0,1,2]);
    tbasisrec:=BasisRec(tGram, 3);
    tOGrec:=NewOGLat(tbasisrec);
    CheckOGrecSize(tbasisrec, tOGrec);
    thesize:=tOGrec.size;
    ttbasisrecs:=[];
    for uuu in [1..5] do 
      tU:=RandomUnimodMat(tdim);
      ttGram:=TMTTmult(tU, tGram);
      ttbasisrec:=BasisRec(ttGram, 3);
      ttOGrec:=NewOGLat(ttbasisrec);
      CheckOGrecSize(ttbasisrec, ttOGrec);
      if ttOGrec.size<>thesize then beep(716211); fi;
      if FindIsomBasisRecs(tbasisrec, ttbasisrec)=fail then beep(919111); fi;
      Add(ttbasisrecs, ttbasisrec);
    od;
    for ss in [1..3] do
      ttbasisrec1:=Random(ttbasisrecs);
      ttbasisrec2:=Random(ttbasisrecs);
      if FindIsomBasisRecs(tbasisrec, ttbasisrec)=fail then beep(229111); fi;
    od;
    Printn(tdim, thesize, tbasisrec.basisnrms);
  od;
od;

["A1", "A2", "A3", "A4","A5", "A6", "A7", "A8", "A9", "A10", 
"D4", "D5", "D6", "D7", "D8", "D9", "D10", "E6", "E7", "E8"];



#####