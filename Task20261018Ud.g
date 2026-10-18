
#Read("Task20261018Ud.g");

dones:=[];
orbs:=[];

for ii in [1..Length(spconsset)] do 
  if not ii in dones then 
    neworb:=Set(t2gperms, tgp->ii^tgp);
    Add(orbs, neworb);
    Append(dones, neworb);
  fi;
od;

Collected(List(orbs, Length));

SortParallel(List(orbs, Length), orbs);
orbs:=Reversed(orbs);

for torb in orbs do 
  Printn(Length(torb), SpconsIntTypes(tGramS, List(torb, ii->spconsset[ii])));
od;


#####

Aut4455rec:=rec(
  AutX:=Aut4455X, 
  spconsset:=spconsset, 
  permsspconsset:=t2gperms,
  orbs:=orbs
);;

savedata(Aut4455rec);



######