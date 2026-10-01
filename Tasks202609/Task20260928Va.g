#Read("Task20260928Va.g");

#
#check NewOGLat and FindIsomBasisRecs
#

#Read("gapReads/OGLat.g");

maxsizes:=List([4..15], ii->0);
for ttt in [1..10000] do
  for tdim in [4..15] do
    tGram:=RandomPosLat(tdim, [-2,-1,0,0,0,0, 1,2]);
    tbasisrec:=BasisRec(tGram, 30);
    tOGrec:=OGLat(tGram);
    CheckOGrecSize(tbasisrec, tOGrec);
    thesize:=tOGrec.size;
    ttbasisrecs:=[];
    for uuu in [1..7] do 
      tU:=RandomUnimodMat(tdim);
      ttGram:=TMTTmult(tU, tGram);
      ttbasisrec:=BasisRec(ttGram, 30);
      ttOGrec:=OGLat(ttGram);
      CheckOGrecSize(ttbasisrec, ttOGrec);
      if ttOGrec.size<>thesize then beep(716211); fi;
      if IsIsomLats(tGram, ttGram)=false then beep(919111); fi;
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