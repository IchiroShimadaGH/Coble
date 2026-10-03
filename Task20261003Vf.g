#Read("Task20261003Vf.g");


readdata("GramLeech");
#Leech4s:=ShortestVectors(GramLeech, 4).vectors;;
#savedata(Leech4s);
#readdata("Leech4s");

resetMT(sessionnumb);

Read("SplitConsTools.g");

RandomFound:=[];
tname:=Concatenation("RandL", String(rho), "H", String(sessionnumb));

sntA:=[];#ANNW


conts:=[0,0,0];

chnops:=24-(rho-1);
tcc:=0;

while true do 
  tcc:=tcc+1;
  if tcc mod 10=0 then Printn(tcc, conts); fi;
  tvs:=RandomChooseFromL(Leech4s, chnops);
  if Rank(tvs)<chnops then conts[1]:=conts[1]+1;  continue; fi;
  tA:=TMTTmult(tvs, GramLeech);
  if DeterminantIntMat(tA)>500 then continue; fi;
  sntA:=SmithNormalFormIntegerMat(tA);
  diags:=List([1..chnops-1], ii->sntA[ii][ii]);
  nzs:=Filtered(diags, xx-> xx<>1);
  if Length(nzs)+1>22-rho then conts[2]:=conts[2]+1;  continue; fi;
  tG:=OrthogonalCompRec(GramLeech, tvs).Gram;
  nopscons:=Length(ShortestVectors(tG, 4).vectors);
  #
  tGramS:=DiagonalMats([ [[2]], -tG]);
  K3flag:=PrimitivelyEmbeddableInK3Lattice(tGramS); 
  if not K3flag then conts[3]:=conts[3]+1;  continue; fi;
  #
  th:=MakeVectei(rho, 1);
  sprats:=GetSpRats (tGramS, th);
  if sprats.nopss[2]<>nopscons then beep(6186168); fi;
  types:=SpconsIntTypes(tGramS, sprats.scons);
  discS:=2*DeterminantIntMat(tG);
  Printn(rho, sprats.nopss, discS, types, tcc);
  Add(RandomFound, [sprats, tvs]) ;
  savedataas(RandomFound, tname);
od;
