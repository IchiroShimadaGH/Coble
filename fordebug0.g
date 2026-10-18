 #Read("Task20261018Ub.g");



tgs:=the4455rec.hperprec.AutShgens;;
tGramS:=the4455rec.GramS;
th:=the4455rec.h;
spcons:=the4455rec.sprats.scons;


involP2:=DiagonalMats([ [[1]], -IdentityMat(19)]);
if th*involP2<>th then beep(332214); fi;
if TMTTmult(involP2, tGramS)<>tGramS then beep(992912); fi;

Add(tgs, involP2);

spconss:=Concatenation(spcons);
Sort(spconss);

if Rank(spconss)<>20 then beep(9919192); fi;

tgperms:=[];
AutX:=[];

for tg in tgs do 
  if th*tg<>th then beep(33214); fi;
  if TMTTmult(tg, tGramS)<>tGramS then beep(99912); fi;
  tgperm:=PermList(List(spconss*tg, tv->PositionSorted(spconss, tv)));
  Add(tgperms, tgperm);
od;

if Size(Group(tgperms))<>11520 then beep(919192); fi;

totalgs:=[IdentityMat(20)];

ii:=0;

while ii<Length(totalgs) do 
  ii:=ii+1;
  tgii:=totalgs[ii];
  for tgen in tgs do 
    newg:=tgii*tgen;
    if not  newg in totalgs then 
      Add(totalgs, newg);
    fi;
  od;
  if ii mod 100=0 then Printn(ii, Length(totalgs)); fi;
od;

if Length(totalgs)<>11520 then beep(988192); fi;


savedata(totalgs);

 ##############