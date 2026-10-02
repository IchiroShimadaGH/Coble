#Read("Task20261003Va.g");


tG:=4*IdentityMat(19);

results:=[];
sntG:=[];#ANNW

while true do 
  ii:=Random([1..19]);
  jj:=Random([1..19]);
  if ii=jj then continue; fi;
  if tG[ii][jj]=0 then 
    aa:=Random([1,-1]);
    tG[ii][jj]:=aa;
    tG[jj][ii]:=aa;
  else 
    aa:=tG[ii][jj];
    tG[ii][jj]:=-aa;
    tG[jj][ii]:=-aa;
  fi;
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
      DeterminantIntMat(tGramS));
      Add(results, [tGramS, th, sprats, List(tG)]) ;
      savedata(results);
    fi;
  fi;
   
od;
