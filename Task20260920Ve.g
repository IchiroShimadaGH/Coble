

#Read("Task20260920Ve.g");

resetMT(sessionnumb);

Read("SplitConsTools.g");

specialadjs2:=[];

TT:=[[-2]];

task:=function(adj, ss)
  local kk, tA, ttA, newadj, trec,tGram, th,isgeom, scons,
  cc, tss;
  kk:=Length(adj);
  for cc in [1..ss] do
    tA:=[RandomVectFromL(kk, [0,1,2,3,4])];
    ttA:=TransposedMat(tA);
    newadj:=MatMatToMat([[adj, ttA], [tA, TT]]);
    trec:=WGraphToGramh(newadj, kk+1);
    tGram:=trec.Gram;
    th:=trec.h;
    isgeom:=IsGeom(tGram, th);
    if isgeom=true then 
      scons:=GetTotalSplcons(tGram, th);
      tss:=Length(scons);
      if tss>120 then 
        Printn(kk+1, Length(tGram), tss);
        Add(specialadjs2, newadj);
        savedata(specialadjs2);
      fi;
      task(newadj, tss);
    fi;
  od;
end;

while true do
 task([[2]], 10);
od;



############