#Read("Task20260909Ud.g");

#
#check NewOGLat and FindIsomBasisRecs
#

Read("NewOGLat.g");

maxsizes:=List([4..15], ii->0);
for ttt in [1..10000] do
  for tdim in [4..15] do
    tGram:=RandomPosLat(tdim, [-2,-1,0,0,0,0, 1,2]);
    tbasisrec:=BasisRec(tGram, 30);
    tOGrec:=NewOGLat(tbasisrec);
    CheckOGrecSize(tbasisrec, tOGrec);
    thesize:=tOGrec.size;
    ttbasisrecs:=[];
    for uuu in [1..5] do 
      tU:=RandomUnimodMat(tdim);
      ttGram:=TMTTmult(tU, tGram);
      ttbasisrec:=BasisRec(ttGram, 30);
      ttOGrec:=NewOGLat(ttbasisrec);
      CheckOGrecSize(ttbasisrec, ttOGrec);
      if ttOGrec.size<>thesize then beep(716211); fi;
      if IsIsomBasisRecs(tbasisrec, ttbasisrec)=false then beep(919111); fi;
      Add(ttbasisrecs, ttbasisrec);
    od;
    for ss in [1..3] do
      ttbasisrec1:=Random(ttbasisrecs);
      ttbasisrec2:=Random(ttbasisrecs);
      if IsIsomBasisRecs(tbasisrec, ttbasisrec)=false then beep(229111); fi;
    od;
    Printn(ttt, tdim, thesize, tbasisrec.basisnrms);
    tdimpos:=SinglePosition([4..15], tdim);
    maxsizes[tdimpos]:=Maximum(maxsizes[tdimpos], thesize);
  od;
  Printn("_____________", maxsizes);
od;




#####