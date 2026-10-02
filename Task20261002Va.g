
#Read("Task20261002Va.g");


#Read("SplitConsTools.g");

tG:=(4)*IdentityMat(19);
for ii in [1..19] do 
  for jj in [1..19] do 
    if AbsInt(ii-jj)=1 then tG[ii][jj]:=1;fi;
  od;
od;

tGramS:=DiagonalMats([ [[2]], -tG]);
th:=MakeVectei(20, 1);

if IsGeom(tGramS, th)<>true then beep(58717861); fi;

sprats:=GetSpRats(tGramS, th);

Printn(sprats.nopss, SpconsIntTypes(tGramS, sprats.scons));
