
 #Read("Task20261016Ua.g");

Read("SplitConsTools.g");

readdata("GramLeech");
#readdata("Leech4s");

kk:=5;

TCC:=NrCombinations([1..24], 5);
iter:=IteratorOfCombinations([1..24], kk);
 
singK3configs:=[];

discgs:=[];
sntA:=[]; #ANNW

tcc:=0;


SpconsIntType:=function(tGram, spcon1, spcon2)
  local tM;
  tM:=spcon1*tGram*TransposedMat(spcon2);
  if tM=[[4,0], [0,4]] or tM=[[0,4], [4,0] ] then return("zf"); fi;
  if tM=[[1,3], [3,1]] or tM=[[3,1], [1,3] ] then return("odd"); fi;
  if tM=[[2,2], [2,2]] then return("tt"); fi;
  beep(1142444111);
end;

FastSpconsIntTypes:=function(tGram, spcons)
  local types, pos1, pos2, spconsdual,  Tspconsdual, nops, spcon1, spcon2,
  tM, ftM, spcon1dual, Tspcon2;
  nops:=Length(spcons);
  Printn("___nops", nops);
  spconsdual:= List(spcons, xx->xx*tGram);
  Tspconsdual:= List(spcons, TransposedMat);
  types:=[0,0,0];
  for pos1 in [1..nops] do
    spcon1dual:=spconsdual[pos1];
    for pos2 in [pos1+1..nops] do
      Tspcon2:=Tspconsdual[pos2];
      tM:=spcon1dual*Tspcon2;
      ftM:=Set(Flat(tM));
      if ftM=[0,4] then types[1]:=types[1]+1; 
      elif ftM=[2] then types[2]:=types[2]+1;  
      elif ftM=[1,3] then types[3]:=types[3]+1; 
      else beep(5875871);
      fi;
    od;
  od;
  return(types);
end;



for tt in iter do 
  tcc:=tcc+1;
  tA:=SubMatrix(GramLeech, tt, tt);
  detA:=DeterminantIntMat(tA);
  if detA>=250 then  continue; fi;
  sntA:=SmithNormalFormIntegerMat(tA);
  diags:=List([1..kk], ii->sntA[ii][ii]);
  notones:=Filtered(diags, xx-> xx<>1);
  if Length(notones)>2 then continue; fi;
  tvs:=List(tt, jj->MakeVectei(24, jj));
  tG:=OrthogonalCompRec(GramLeech, tvs).Gram;
  th:=MakeVectei(20, 1);
  tGramS:=DiagonalMats([ [[2]], -tG]);
  discg:=DiscriminantForm(tGramS).discg;
  if discg in discgs then continue; fi;
  Add(discgs, discg);
  #
  K3flag:=PrimitivelyEmbeddableInK3Lattice(tGramS); 
  if not K3flag then 
    Printn(__tt);
    continue; 
  fi;
  nopscons:=Length(ShortestVectors(tG, 4).vectors);
  sprats:=GetSpRats(tGramS, th);
  if sprats.nopss[2]<>nopscons then beep(74111); fi;
  dd:=AbsInt(DeterminantIntMat(tGramS));
  if dd<>2*detA then beep(55111); fi;
  if sprats.slins<>[] then beep(616818); fi;
  gtypes:=FastSpconsIntTypes(tGramS, sprats.scons);
  #
  Printn(tt, dd, sprats.nopss, gtypes);
  #
  trec:=rec(
    tt:=tt, 
    rr:=Length(tGramS), 
    dd:=dd, 
    GammaSize:=sprats.nopss[2], 
    LambdaSize:=sprats.nopss[1], 
    gtypes:=gtypes,
    GramS:=tGramS,
    h:=th,
    sprats:=sprats
  );
  #
  Add(singK3configs, trec);
od;

dds:=List(singK3configs, xx->xx.dd);
SortParallel(dds, singK3configs);

for xx in singK3configs do 
  Printn(xx.dd, 1, xx.GammaSize, xx.LambdaSize, xx.gtypes);
od;

savedata(singK3configs);


#######################