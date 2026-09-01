

aR:=2*Sum(GramRdual);

orthorec:=OrthogonalCompRec(GramR, [aR]);
if ShortestVectors(-orthorec.Gram, -2).vectors<>[] then beep(112211); fi;


aRL:=aR*embR;
aL:=aRL;

flag:=false;
while not flag  do 
  aL:=aL+aSL;
  flag:=() and (AffESstd(GramL, aL, 0, -2, true)=[])
od;

if AffESstd(GramL, aL, 0, -2, true)<>[] then beep(42231); fi;

 Lroots:=[];

dd:=0;
 while Lroots=[] or Rank(Lroots)<26 do 
  dd:=dd+1;
  newrs:=Add(AffESstd(GramL, aL, dd, -2, true));
  Append(Lroots, newrs);
  Printn(dd, Length(newrs), Rank(Lroots));
 od;

aL*GramL*aL;