
#Read("Task20261001Ub.g");

Read("NonSingOverLats.g");
Read("even_lattice_genus.g");
Read("SplitConsTools.g");

readdata("Data1845");



Gram:=Data1845.Gram;
 
spcons:=Data1845.scons;

if Length(spcons)<>1845 then beep(581815); fi;

types:=[];

for pos1 in [1..1845] do
  spcon1:=spcons[pos1];
  for pos2 in [pos1+1..1845] do
    spcon2:=spcons[pos2];
    Add(types, SpconsIntType(Gram, spcon1, spcon2));
  od;
  Printn(pos1, Collected(types));
od;

Printn(Collected(types));


##########

