#Read("Task20260909Ua.g");


#
#check NewOGLat and FindIsomBasisRecs
#

Read("NewOGLat.g");

difficult12:=[];

for ttt in [1..10000] do
  tGram3:=RandomPosLat(3, [-4,-3, -2,-1,0,0,0,0, 1,2, 3,4]);
  tbasisrec3:=BasisRec(tGram3, 30);
  tOGrec3:=RepeatNewOGLat(tGram3, 10, 10);
  tsize3:=tOGrec3.size;
  themm:=24*tsize3^4;
  tGram3s:=[];
  for kk in [1..4] do 
    tU:=RandomUnimodMat(3);
    Add(tGram3s, TMTTmult(tU, tGram3));
  od;
  tGram12:=DiagonalMats(tGram3s);
  #
  st:=Runtime();
  tU:=RandomUnimodMat(12);
  tbasisrec:=BasisRec(tGram12, 30);
  tOGrec:=RepeatNewOGLat(tGram12, 10, 10);
  CheckOGrecSize(tbasisrec, tOGrec);
  thesize:=tOGrec.size;
  if thesize mod themm <>0 then beep(888111); fi;
  ttime:=Runtime()-st;
  if ttime>1000*100 then 
    Add(difficult12, tGram12);
    savedata(difficult12);
  fi;
  Printn("____", thesize, themm, tbasisrec.basisnrms, TimeToString(ttime));
  ttbasisrecs:=[];
  for uuu in [1..5] do 
    tU:=RandomUnimodMat(12);
    ttGram:=TMTTmult(tU, tGram12);
    ttbasisrec:=BasisRec(ttGram, 30);
    ttOGrec:=RepeatNewOGLat(ttGram, 10, 10);
    CheckOGrecSize(ttbasisrec, ttOGrec);
    if ttOGrec.size<>thesize then beep(712211); fi;
    fflag:=IsIsomBasisRecs(tbasisrec, ttbasisrec);
    if fflag=false then beep(919111); fi;
    if not IsIntMat(fflag) then beep(919331); fi;
    Add(ttbasisrecs, ttbasisrec);
    Printn("______", uuu, ttOGrec.size, ttbasisrec.basisnrms);
  od;
  for ss in [1..3] do
    ttbasisrec1:=Random(ttbasisrecs);
    ttbasisrec2:=Random(ttbasisrecs);
    fflag:=IsIsomBasisRecs(tbasisrec, ttbasisrec);
    if fflag=false then beep(922111); fi;
    if not IsIntMat(fflag) then beep(955332); fi;
  od;
  Printn(ttt, tsize3, thesize, tbasisrec.basisnrms);
  Printn("_____________");
od;




#####





