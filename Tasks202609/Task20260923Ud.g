#Read("Task20260923Ud.g");

Read("SplitConsTools.g");
Read("even_lattice_genus.g");

readdata("Data1845");

adj1845:=Data1845.adj;

kk:=Length(adj1845);  #=17

TT:=[[-2]];

Over1845NotOdd:=[];

counter:=0;
while true do
  counter:=counter+1;
  tA:=[RandomVectFromL(kk, [0,1,2,3,4])];
  ttA:=TransposedMat(tA);
  newadj:=MatMatToMat([[adj1845, ttA], [tA, TT]]);
  trec:=WGraphToGramh(newadj, kk+1);
  tGram:=trec.Gram;
  th:=trec.h;
  isgeom:=IsGeom(tGram, th);
  if isgeom=true then
    scons:=GetTotalSplcons(tGram, th);
    if Length(scons)<>nopsscons then beep(9192992); fi;
    ttmdiscrec:=DiscriminantForm(-tGram);
    tsign:=[22, 3, 19]-SignatureQ(tGram);
    ttsign:=[tsign[2], tsign[3]];
    Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
    if Trec.count=0 then buzz(47652);fi;
    trec.scons:=scons;
    trec.Trec:=Trec;
    if Length(scons)<1845 then buzz(118181); fi;
    Add(Over1845NotOdd, trec);
    Printn(counter, ":", Length(tGram), Length(scons), ttmdiscrec.discg);
    savedata(Over1845NotOdd);
  fi;
  if counter mod 10000 =0 then Printn("___", counter); fi;
od;









############