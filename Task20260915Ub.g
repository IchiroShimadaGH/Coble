#Read("Task20260915Ub.g");

tGram:=[
  [2,2,2,2,2],
  [2, -2,3,1,3],
  [2,3,-2,3,1],
  [2,1,3,-2,3],
  [2,3,1,3,-2]
];

if PrimitivelyEmbeddableInK3Lattice(tGram)<>true then buzz(665511); fi;

th:=[1,0,0,0,0];

if AffESstd(tGram, th, 0, -2, true)<>[] then beuzz(77152); fi;

rats2:=FindRatsOnK3(tGram, th,2);

Sort(rats2);

tg:=RatInvolRecP2(tGram, th, th).invol;

pos:=0;
for tv in rats2 do 
  pos:=pos+1;
  tvtg:=tv*tg;
  if tvtg*tGram*tv<>6 then beep(66161); fi;
  Printn(tv, tv=tvtg, SinglePosition(rats2, tv*tg));
od;
#####