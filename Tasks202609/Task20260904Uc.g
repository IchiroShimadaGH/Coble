#Read("Task20260904Uc.g");

Read("GeneralStabChain.g");
Read("Make_nSaSvSss.g");
Read("Adj.g");

readdata("Borrec");
readdata("nSaSvSss0");
readdata("adjrecs0");

#autdiscRrec:=AutDiscfByGeneralStabChain(Borrec.discR);;
Printn("autdiscRrec.size", autdiscRrec.size);
savedatap(autdiscRrec);

uvabdiscR:=uvabvect(Borrec.discR);
uvabdiscS:=uvabvect(Borrec.discS);
#tsizeR:=uvabOqLOrder(uvabdiscR);
#tsizeS:=uvabOqLOrder(uvabdiscS);
if tsizeR<>autdiscRrec.size then beep(88181); fi;
if tsizeS<>autdiscRrec.size then beep(88181); fi;

OGRrec:=OGLat(-Borrec.GramR);;
savedatap(OGRrec);
if OGRrec.order<> WeylGroupOrder(Borrec.Rtype)*Factorial(7)*2*6*6 then beep(91919); fi;

discgR:=Borrec.discR.discg;

vs:=Cartesian(List(discgR, kk->[0..kk-1]));;
tqgperms:=[];

for tg in OGRrec.totalgens do 
  tqg:=OLtoOqL(tg, Borrec.discR);
  tqgperm:=OqLtoOqLPerm(tqg, discgR, vs);
  Add(tqgperms, tqgperm);
od;

if Size(Group(tqgperms))<>Factorial(7)*2*6*6 then beep(199191); fi;


####################