#Read("Task20261003Ua.g");

resetMT(sessionnumb);

Read("SplitConsTools.g");

randinitG:=function()
  local tG, ii, jj, aa;
  tG:=4*IdentityMat(rho-1);
  for ii in [1..rho-1] do 
    for jj in [ii+1..rho-1] do 
      aa:=Random([-1, 0,0, 0, 1]);
      tG[ii][jj]:=aa;
      tG[jj][ii]:=aa;
    od;
  od;
  return(tG);
end;



RandomFound:=[];
tname:=Concatenation("Rand", String(rho), "Pro", String(sessionnumb));

sntG:=[];#ANNW

tcc:=0;
fcc:=0;
abs:=[0,0,0,0,0,0];
tG:=4*IdentityMat(rho-1);

while true do 
  tcc:=tcc+1;
  fcc:=fcc+1;
  if fcc=20000 then 
    Printn("reset", tcc, "rho", rho, abs);
    fcc:=0;
    tG:=4*IdentityMat(rho-1);
  fi;
  ii:=Random([1..rho-1]);
  jj:=Random([1..rho-1]);
  if ii=jj then abs[1]:=abs[1]+1;  continue; fi;
  aa:=Random([-1..1]);
  tG[ii][jj]:=aa;
  tG[jj][ii]:=aa;
  #
  if Rank(tG)<rho-1 then abs[2]:=abs[2]+1;  continue; fi;
  if SignatureQ(tG)<>[rho-1,rho-1, 0] then abs[3]:=abs[3]+1;  continue; fi;
  sntG:=SmithNormalFormIntegerMat(tG);
  diags:=List([1..rho-1], ii->sntG[ii][ii]);
  nzs:=Filtered(diags, xx-> xx<>1);
  if Length(nzs)>21-rho then abs[4]:=abs[4]+1;  continue; fi;
  svs:=ShortestVectors(tG, 4);
  if 2 in svs.norms then abs[5]:=abs[5]+1;  continue; fi;
  nopscons:=Length(svs.norms);
  #
  tGramS:=DiagonalMats([ [[2]], -tG]);
  K3flag:=PrimitivelyEmbeddableInK3Lattice(tGramS); 
  if not K3flag then abs[6]:=abs[6]+1;  continue; fi;
  #
  th:=MakeVectei(rho, 1);
  sprats:=GetSpRats (tGramS, th);
  if sprats.nopss[2]<>nopscons then beep(6186168); fi;
  types:=SpconsIntTypes(tGramS, sprats.scons);
  discS:=2*DeterminantIntMat(tG);
  if nopscons>=rho then 
    Printn(sprats.nopss, discS, types, tcc, fcc);
    fcc:=0;
    Add(RandomFound, [sprats, List(tG)]) ;
    savedataas(RandomFound, tname);
  fi;
od;
