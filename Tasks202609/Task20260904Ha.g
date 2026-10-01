#Read("Task20260904Ha.g");


Stype:=["H2"];
for ii in [1..10] do 
  Add(Stype, "A1");
od;

discS:=NegADEHdiscrec(Stype);

diagsS:=[2];
for ii in [1..10] do 
  Add(Stype, -2);
od;

GramS:=DiagonalMats(diagsS)

Rtype:=["D6"];
for ii in [1..9] do 
  Add(Rtype, "A1");
od;

discR:=NegADEHdiscrec(Rtype);
mdiscR:=PosADEHdiscrec(Rtype);

GramD6:=ADEHGram(["D6"]);
discD6:=DiscriminantForm(GramD6);
uvabvect(discD6);


a4:=[[2,2,2,2], (1/2)*IdentityMat(4)];
b4:=[[2,2,2,2], (3/2)*IdentityMat(4)];
anisom:=AnIsomFQF(a4, b4);
autb4:=AutFQF(b4);




#####