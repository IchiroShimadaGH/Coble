#Read("Task20261004Va.g");

Read("SplitConsTools.g");

readdata("GramLeech");
#Leech4s:=ShortestVectors(GramLeech, 4).vectors;;
#savedata(Leech4s);
#readdata("Leech4s");


kk:=6;
TCC:=NrCombinations([1..24], 6);

iter:=IteratorOfCombinations([1..24], kk);
 
rho19results:=[];

discgs:=[];

tcc:=0;
for tt in iter do 
  tcc:=tcc+1;
  tA:=SubMatrix(GramLeech, tt, tt);
  detA:=DeterminantIntMat(tA);
  if detA>=250 then  continue; fi;
  sntA:=SmithNormalFormIntegerMat(tA);
  diags:=List([1..kk], ii->sntA[ii][ii]);
  notones:=Filtered(diags, xx-> xx<>1);
  if Length(notones)>3 then continue; fi;
  tvs:=List(tt, jj->MakeVectei(24, jj));
  tG:=OrthogonalCompRec(GramLeech, tvs).Gram;
  tGramS:=DiagonalMats([ [[2]], -tG]);
  discg:=DiscriminantForm(tGramS).discg;
  if discg in discgs then continue; fi;
  Add(discgs, discg);
  #
  K3flag:=PrimitivelyEmbeddableInK3Lattice(tGramS); 
  if not K3flag then  continue; fi;
  nopscons:=Length(ShortestVectors(tG, 4).vectors);
  
  #
  Printn(tt, detA, nopscons, diags, discg,  K3flag, tcc, TCC);
  Add(rho19results, tt);
od;

savedata(rho19results);



