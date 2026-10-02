

#Read("Task20261001Pa.g");

#Random Odd Enhanced 2

Read("SplitConsTools.g");

resetMT(sessionnumb);

specialadjsf:=[];
savename:=Concatenation("randomPro", String(sessionnumb));

TT:=[[-2]];

taskcounter:=0;

task:=function(adj)
  local kk, tA, ttA, newadj, trec,tGram, th,isgeom, srats,
  cc, mm;
  taskcounter:=taskcounter+1;
  kk:=Length(adj);
  mm:=Maximum(50, 50*kk);
  for cc in [1..mm] do
    tA:=[RandomVectFromL(kk, [1, 3])];
    ttA:=TransposedMat(tA);
    newadj:=MatMatToMat([[adj, ttA], [tA, TT]]);
    trec:=WGraphToGramh(newadj, kk+1);
    tGram:=trec.Gram;
    th:=trec.h;
    isgeom:=IsGeom(tGram, th);
    if isgeom=true then 
      srats:=GetSpRats(tGram, th);
      if srats.nopss[2]>1000 or Length(tGram)>17 then 
        Printn(kk+1, Length(tGram), srats.nopss, Collected(Flat(newadj)), taskcounter);
        Add(specialadjsf, [srats.nopss, newadj]);
        savedataas(specialadjsf, savename);
      fi;
      task(newadj);
    fi;
  od;
end;

task([[-2]]);



############