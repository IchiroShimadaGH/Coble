

#Read("Task20261018Uc.g");

Aut4455X:=[];

tdiscSXrec:=DiscriminantForm(tGramS);

if tdiscSXrec.discg<>[6, 60] then beep(19198); fi;

periodqgs:=[[[1,0], [0,1]],  [[5, 0], [0, 59]]];

for tg in totalgs do
  tqg:=OLtoOqL(tg, tdiscSXrec);
  if tqg in periodqgs then 
    Add(Aut4455X, tg);
  fi;
od;

Printn(Length(Aut4455X));

if not involP2 in Aut4455X then beep(514231); fi;

for tg in Aut4455X do if tg*involP2<>involP2*tg then beep(771621); fi; od;

spconsset:=Set(spcons, Set);;

t2gperms:=[];

for tg in Aut4455X do 
  tgp:=PermList(List(spconsset, tvs->PositionSorted(spconsset, Set(tvs*tg))));
  Add(t2gperms, tgp);
od;;

savedata(Aut4455X);




# t2gperms:=[];

# for tg in Aut4455X do 
#   tgperm:=PermList(List(spconss*tg, tv->PositionSorted(spconss, tv)));
#   Add(t2gperms, tgperm);
# od;
# gap> StructureDescription(Group(t2gperms));
# "C2 x A6"


####