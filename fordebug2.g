
K16gens:=erecs[3].minusrec.Kerqrec.gens;;

basism:=erecs[3].minusrec.basis;
Gramm:=erecs[3].minusrec.Gram;

tau120:=SolutionIntMat(basism, gammarecs[4].tau);
tau135:=SolutionIntMat(basism, gammarecs[5].tau);

Xitau120:=List(gammarecs[4].plusrec.Xitau, tv->SolutionMat(basism, tv));
Xitau135:=List(gammarecs[5].plusrec.Xitau, tv->SolutionMat(basism, tv));

tilXitau120:=List(Xitau120, tv->tv-tau120/2);
tilXitau135:=List(Xitau135, tv->tv-tau135/2);

K16permgens120:=[];
K16permgens135:=[];

for tg in K16gens do 
  tgperm:=PermList(List(tilXitau120*tg, tv->SinglePosition(tilXitau120, tv))); 
  Add(K16permgens120, tgperm);
  tgperm:=PermList(List(tilXitau135*tg, tv->SinglePosition(tilXitau135, tv))); 
  Add(K16permgens135, tgperm);
od;


if Size(Group(K16permgens120))<>512 then beep(73342); fi;
if Size(Group(K16permgens135))<>512 then beep(776542); fi;

if List(OrbitdDecompByPerms(Length(tilXitau120), K16permgens120), Length)<>[512] then 
  beep(58185185); 
fi;
if List(OrbitdDecompByPerms(Length(tilXitau135), K16permgens135), Length)<[32] then 
  beep(58185185); 
fi;

for tt in [1..50] do
  aa:=Random([1..512]);
  if Set(Group(K16permgens120), tp->aa^tp)<>[1..512] then beep(88811); fi;
od;


for aa in [1..32] do
  if Set(Group(K16permgens135), tp->aa^tp)<>[1..32] then beep(22811); fi;
od;


#############