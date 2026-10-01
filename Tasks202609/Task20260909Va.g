#Read("Task20260909Va.g");


readdata("Rrecs");

for trec in Rrecs do 
  tGram:=trec.Gram;
  for tg in trec.gens do 
    if TMTTmult(tg, tGram)<>tGram then beep(999111); fi;
  od;
od;


GensOGSize:=function(tGram, gens)
  local vs, tg, tperm, gensperms, permsize;
  vs:=ShortestVectors(tGram, 4).vectors;
  vs:=Union(vs, -vs);
  gensperms:=[];
  for tg in gens do
    tperm:=PermList(List(vs*tg, ttv->Position(vs, ttv)));
    Add(gensperms, tperm);
  od;
  permsize:=Size(Group(gensperms));
  return(permsize);
end;

pos:=0;
for trec in Rrecs do 
  pos:=pos+1;
  if trec.ogsize<>GensOGSize(trec.Gram, trec.gens) then 
    beep(81818);
  else 
    Printn("ogsize check", pos);
  fi;
od;


pos:=0;
for trec in Rrecs do 
  pos:=pos+1;
  imsize:=ImOLtoOqL(trec.gens, trec.Gram);
  kersize:=KerOLtoOqL(trec.Gram).size;
  Unbind(trec.imetaindex);
  trec.kersize:=kersize;
  trec.imsize:=imsize;
  if trec.ogsize<>imsize*kersize then 
    beep(8111818);
  else 
    Printn("imetasize check", pos, trec.imsize, Collected(trec.ADEtype));
  fi;
od;

savedata(Rrecs);

####################