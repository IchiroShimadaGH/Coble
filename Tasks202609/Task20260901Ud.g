#Read("Task20260901Ud.g");

readdata("SRLdata");
readdata("LampleS");

GramS:=SRLdata.GramS;
GramR:=SRLdata.GramR;
GramL:=SRLdata.GramL;
embS:=SRLdata.embS;
embR:=SRLdata.embR;
projS:=SRLdata.projS;
projR:=SRLdata.projR;


aS:=LampleS;
hS:=MakeVectei(11,1);
if SeparatingVects(GramS, hS, aS, -2)<>[] then beep(81811); fi;

aR:=2*Sum(InverseMat(GramR));

orthorec:=OrthogonalCompRec(GramR, [aR]);
if ShortestVectors(-orthorec.Gram, -2).vectors<>[] then beep(112211); fi;

rootsR:=ShortestVectors(-GramR, 2).vectors;
rootsR:=Union(rootsR, -rootsR);

if not IsEqualSet(AffESstd(GramL, aSL, 0, -2, true), rootsR*embR) then 
  beep(991919); 
fi;

aSL:=aS*embS;
aRL:=aR*embR;
aL:=aRL;

flag:=false;
while not flag  do 
  aL:=aL+aSL;
  dd:=aL*GramL*aL;
  if dd>0 then 
    tvs:=AffESstd(GramL, aL, 0, -2, true);
    Printn(dd, Length(tvs));
    if tvs=[] then 
      flag:=true; 
      seps:=SeparatingVects(GramL, aSL, aL, -2);
      if seps<>[] then 
        minti:=seps[1][2];
        ttv:=aSL+(minti/2)*aL;
        mm:=Llcm(List(ttv, DenominatorRat));
        if mm<=0 then beep(5875187); fi;
        aL:=mm*ttv;
      fi;
    fi;
 fi;
od;

if not IsIntVect(aL) then beep(6986198); fi;
if AffESstd(GramL, aL, 0, -2, true)<>[] then beep(229819); fi;
if SeparatingVects(GramL, aSL, aL, -2)<>[] then beep(919819); fi;


amplesdata:=rec(
  hS:=hS,
  hSL:=hS*embS, 
  aS:=aS,
  aSL:=aSL,
  aR:=aR, 
  aRL:=aRL,
  aL:=aL
);

savedata(amplesdata);


##################