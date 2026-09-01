#Read("Task20260901Ub.g");

aa:=List([1..10], ii->-2);
Add(aa, 2, 1);

GramS:=DiagonalMat(aa);
discS:=DiscriminantForm(GramS);
discmS:=DiscriminantForm(-GramS);

Rtype:=Concatenation(["A1", "A1", "A1", "A1", "A1", "A1", "A1"], ["D4", "D4"]);
GramR:=-ADEHGram(Rtype);
discR:=DiscriminantForm(GramR);

# anisom:=AnIsomFQF([discmS.discg, discmS.discf], [discR.discg, discR.discf]);

# Printn(anisom);

# savedata(anisom);
# savedataas(anisom, "anisomfromdiscmStodiscR");

#readdata("anisom");

if CopyNormalDiscf(TMTTmult(anisom, discR.discf))<>discmS.discf then beep(686827); fi;


GramSdual:=InverseMat(GramS);
GramRdual:=InverseMat(GramR);
GramLbonedual:=DiagonalMats([GramSdual, GramRdual]);
tgraph:=ConcatMats(discS.reps_dual, anisom*discR.reps_dual);
if not IsIntMat(TMTTmult(tgraph, GramLbonedual)) then beep(996918); fi;

bonebasis:=DiagonalMats([GramS, GramR]);
Lgens:=StackMats(bonebasis,tgraph);
hLgens:=HermiteNormalFormIntegerMat(Lgens);
if not IsZeroMat(List([27..Length(hLgens)], ii->hLgens[ii])) then beep(721711); fi;
tLbasis:=List([1..26], ii->hLgens[ii]);

GramL:=TMTTmult(tLbasis, GramLbonedual);
if not IsEvenLattice(GramL) then beep(661961); fi;
if SignatureQ(GramL)<>[26,1,25] then beep(12345); fi;
if DeterminantIntMat(GramL)<>-1 then beep(919191); fi;

embS:=List(ConcatMats(GramS, NullMat(11, 15)), tv->SolutionIntMat(tLbasis, tv));
embR:=List(ConcatMats(NullMat(15, 11), GramR), tv->SolutionIntMat(tLbasis, tv));

if TMTTmult(embS, GramL)<>GramS then beep(9998111); fi;
if TMTTmult(embR, GramL)<>GramR then beep(9333111); fi;

if CokerTorsion(embS)<>[] then beep(888111); fi;
if CokerTorsion(embR)<>[] then beep(822211); fi;

embsinv:=InverseMat(StackMats(embS, embR));

projS:=SubMatrix(embsinv, [1..26], [1..11]);
projR:=SubMatrix(embsinv, [1..26], [12..26]);

if embS*projS<>Idmat11 then beep(919191); fi;
if embR*projR<>Idmat15 then beep(229191); fi;
if not IsZeroMat(embS*projR) then beep(56198); fi;
if not IsZeroMat(embR*projS) then beep(444498); fi;



hS:=MakeVectei(11, 1);
if hS*GramS*hS<>2 then beep(95911); fi;

aS:=List([1..11], ii->0);
aS[1]:=4;
for kk in [2..11] do 
  aS[kk]:=-1;
od;


if aS*GramS*aS<=0 then beep(44911); fi;
if hS*GramS*aS<=0 then beep(25911); fi;
if AffESstd(GramS, aS, 0, -2, true)<>[] then beep(44431); fi;
if SeparatingVects(GramS, hS, aS, -2)<>[] then beep(22911); fi;

nrR:=ADEHrootnumb(Rtype);

SRLdata:=rec(
  GramS:=GramS, 
  GramR:=GramR, 
  GramL:=GramL, 
  embS:=embS, 
  embR:=embR, 
  projR:=projR, 
  projS:=projS,
  Rtype:=Rtype, 
  discS:=discS,
  discR:=discR,
  anisomdiscSR:=anisom 
);

savedata(SRLdata);


#################