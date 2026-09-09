#Read("Task20260909Ub.g");

#
#check NewOGLat and FindIsomBasisRecs
#

Read("NewOGLat.g");

difficult16:=[];

for ttt in [1..10000] do
  tGram4:=RandomPosLat(4, [-2,-1,0,0,0,0, 1,2]);
  tbasisrec4:=BasisRec(tGram4, 30);
  tOGrec4:=RepeatNewOGLat(tGram4, 10, 10);
  tsize4:=tOGrec4.size;
  themm:=24*tsize4^4;
  tGram4s:=[];
  for kk in [1..4] do 
    tU:=RandomUnimodMat(4);
    Add(tGram4s, TMTTmult(tU, tGram4));
  od;
  tGram16:=DiagonalMats(tGram4s);
  #
  st:=Runtime();
  tbasisrec:=BasisRec(tGram16, 30);
  tOGrec:=RepeatNewOGLat(tGram16, 10, 10);
  CheckOGrecSize(tbasisrec, tOGrec);
  thesize:=tOGrec.size;
  if thesize mod themm <>0 then beep(888111); fi;
  ttime:=Runtime()-st;
  if ttime>1000*100 then 
    Add(difficult16, tGram16);
    savedata(difficult16);
  fi;
  Printn("____", thesize, themm, tbasisrec.basisnrms, TimeToString(ttime));
  ttbasisrecs:=[];
  for uuu in [1..5] do 
    tU:=RandomUnimodMat(16);
    ttGram:=TMTTmult(tU, tGram16);
    ttbasisrec:=BasisRec(ttGram, 30);
    ttOGrec:=RepeatNewOGLat(ttGram, 10, 10);
    CheckOGrecSize(ttbasisrec, ttOGrec);
    if ttOGrec.size<>thesize then beep(716211); fi;
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
  Printn(ttt, tsize4, thesize, tbasisrec.basisnrms);
  Printn("_____________");
od;




#####