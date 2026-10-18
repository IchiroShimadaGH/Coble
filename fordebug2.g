
dones:=[];
orbs:=[];

for ii in [1..Length(spconsset)] do 
  if not ii in dones then 
    neworb:=Set(t2gperms, tgp->ii^tgp);
    Add(orbs, neworb);
    Append(dones, neworb);
  fi;
od;

Collect(List(orbs, Length));
