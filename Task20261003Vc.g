#Read("Task20261003Vc.g");

Read("SplitConsTools.g");
tG:=4*IdentityMat(19);

results2:=[];
sntG:=[];#ANNW

tcc:=0;
fcc:=0;
while true do 
  tcc:=tcc+1;
  fcc:=fcc+1;
  if fcc=10000 then 
    Printn("reset", tcc);
    fcc:=0;
    tG:=4*IdentityMat(19);
  fi;
  ii:=Random([1..19]);
  jj:=Random([1..19]);
  if ii=jj then continue; fi;
  aa:=Random([1, -1]);
  tG[ii][jj]:=aa;
  tG[jj][ii]:=aa;
  if Rank(tG)<19 then continue; fi;
  if SignatureQ(tG)<[19,19, 0] then continue; fi;
  sntG:=SmithNormalFormIntegerMat(tG);
  diags:=List([1..19], ii->sntG[ii][ii]);
  nzs:=Filtered(diags, xx-> xx<>1);
  if Length(nzs)=1 then 
    tGramS:=DiagonalMats([ [[2]], -tG]);
    th:=MakeVectei(20, 1);
    if IsGeom(tGramS, th)=true then
      sprats:=GetSpRats (tGramS, th);
      Printn(sprats.nopss, SpconsIntTypes(tGramS, sprats.scons), 
      DeterminantIntMat(tGramS), tcc, fcc);
      fcc:=0;
      Add(results2, [tGramS, th, sprats, List(tG)]) ;
      savedata(results2);
    fi;
  fi;
   
od;
