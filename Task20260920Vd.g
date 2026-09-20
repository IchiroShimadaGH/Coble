

#Read("Task20260920Vd.g");


Read("SplitConsTools.g");

resetMT(sessionnumb);

specialadjs:=[];

TT:=[[-2]];

task:=function(adj)
  local kk, tA, ttA, newadj, trec,tGram, th,isgeom, scons,
  cc;
  kk:=Length(adj);
  for cc in [1..100] do
    tA:=[RandomVectFromL(kk, [0,1,2,3,4])];
    ttA:=TransposedMat(tA);
    newadj:=MatMatToMat([[adj, ttA], [tA, TT]]);
    trec:=WGraphToGramh(newadj, kk+1);
    tGram:=trec.Gram;
    th:=trec.h;
    isgeom:=IsGeom(tGram, th);
    if isgeom=true then 
      scons:=GetTotalSplcons(tGram, th);
      if Length(scons)>120 then 
        Printn(kk, Length(scons));
        Add(specialadjs, newadj);
        savedata(specialadjs);
      fi;
      task(newadj);
    fi;
  od;
end;

task([[2]]);



############