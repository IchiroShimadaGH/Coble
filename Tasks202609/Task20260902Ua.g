#Read("Task20260902Ua.g");

readdata("SRLdata");
readdata("LampleS");
readdata("amplesdata");

GramS:=SRLdata.GramS;
GramR:=SRLdata.GramR;
GramL:=SRLdata.GramL;
embS:=SRLdata.embS;
embR:=SRLdata.embR;
projS:=SRLdata.projS;
projR:=SRLdata.projR;

hS:=amplesdata.hS;
hSL:=amplesdata.hSL;
aS:=amplesdata.aS;
aSL:=amplesdata.aSL;
aR:=amplesdata.aR;
aRL:=amplesdata.aRL;
aL:=amplesdata.aL;

if hSL<>hS*embS then beep(167161); fi;
if aSL<>aS*embS then beep(267161); fi;
if aRL<>aR*embR then beep(367161); fi;


if AffESstd(GramS, aS, 0, -2, true)<>[]  then beep(382819); fi;
if AffESstd(GramL, aL, 0, -2, true)<>[]  then beep(481819); fi;

if SeparatingVects(GramS, hS, aS, -2)<>[] then beep(191991); fi;
if SeparatingVects(GramL, aSL, aL, -2)<>[] then beep(291991); fi;

ratsL:=[];

dd:=0;
while ratsL=[] or Rank(ratsL)<26 do 
  dd:=dd+1;
  newrats:=AffESstd(GramL, aL, dd, -2, true);
  Append(ratsL, newrats);
  if ratsL<>[] then 
    Printn(dd, Length(newrats), Rank(ratsL));
  else 
    Printn(dd, Length(newrats));
  fi;
od;

savedata(ratsL);

ratsLdualsT:=TransposedMat(ratsL*GramL);

ones:=List([1..Length(ratsL)], ii->1);

weyl0:=SolutionIntMat(ratsLdualsT, ones);

if weyl0=fail then beep(88811); fi;
if weyl0*GramL*weyl0<>0 then beep(991233); fi;
if weyl0*GramL*aL<=0 then beep(66152); fi;

savedata(weyl0);



##################