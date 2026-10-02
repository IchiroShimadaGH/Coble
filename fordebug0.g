#read("Task20261003Va.g");
)

tG:=4*IdentityMat(19);

resuts:=[];

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
  sntG:=SmithNormalFormIntegerMat(tG);
  nzs:=Filtered(List([1..19], ii->sntG[ii][ii], xx-> xx<>1));
  if Length(nzs)=1 then 
    tGramS:=DiagonalMats([ [[2]], -tG]);
    th:=MakeVectei(20, 1);
    if IsGeom(tGramS, th)=true then 
      Add(resuts, List(tG)) ;
    fi;
  fi;
   
od;
