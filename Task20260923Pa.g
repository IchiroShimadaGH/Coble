

#Read("Task20260923Pa.g");

#Random Odd Enhanced 2

Read("SplitConsTools.g");

resetMT(sessionnumb);

specialadjsf:=[];
savename:=Concatenation("randomPro", String(sessionnumb));

TT:=[[-2]];

task:=function(adj)
  local kk, tA, ttA, newadj, trec,tGram, th,isgeom, scons,
  cc, mm;
  kk:=Length(adj);
  mm:=Maximum(50, 30*kk);
  if kk>13 then mm:=3^(kk-13)*mm;fi;
  for cc in [1..mm] do
    tA:=[RandomVectFromL(kk, [1, 3])];
    ttA:=TransposedMat(tA);
    newadj:=MatMatToMat([[adj, ttA], [tA, TT]]);
    trec:=WGraphToGramh(newadj, kk+1);
    tGram:=trec.Gram;
    th:=trec.h;
    isgeom:=IsGeom(tGram, th);
    if isgeom=true then 
      scons:=GetTotalSplcons(tGram, th);
      if Length(scons)>300 then 
        Printn(kk+1, Length(tGram), Length(scons), Collected(Flat(newadj)));
        Add(specialadjsf, [Length(scons), newadj]);
        savedataas(specialadjsf, savename);
      fi;
      task(newadj);
    fi;
  od;
end;

task([[-2]]);



############