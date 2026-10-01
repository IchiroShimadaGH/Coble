

#Read("Task20260922Vb.g");

#Random Even 

Read("SplitConsTools.g");

resetMT(sessionnumb);

specialadjsf:=[];
savename:="randomeven";

TT:=[[-2]];

task:=function(adj)
  local kk, tA, ttA, newadj, trec,tGram, th,isgeom, scons,
  cc;
  kk:=Length(adj);
  for cc in [1..100] do
    tA:=[RandomVectFromL(kk, [0,2,4])];
    ttA:=TransposedMat(tA);
    newadj:=MatMatToMat([[adj, ttA], [tA, TT]]);
    trec:=WGraphToGramh(newadj, kk+1);
    tGram:=trec.Gram;
    th:=trec.h;
    isgeom:=IsGeom(tGram, th);
    if isgeom=true then 
      scons:=GetTotalSplcons(tGram, th);
      if Length(scons)>120 then 
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