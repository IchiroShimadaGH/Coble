#Read("Task20260923Ue.g");

Read("SplitConsTools.g");
Read("even_lattice_genus.g");

Printn("All-zero adj:");
AllZeroAdjRecs:=[];

for kk in [2..19] do 
  adjmat:=(-2)*IdentityMat(kk);
  trec:=WGraphToGramh(adjmat, kk);
  tGram:=trec.Gram;
  th:=trec.h;
  isgeom:=IsGeom(tGram, th);
  trec.isgemo:=isgeom;
  ttmdiscrec:=DiscriminantForm(-tGram);
  tsign:=[22, 3, 19]-SignatureQ(tGram);
  ttsign:=[tsign[2], tsign[3]];
  if tsign[3]>0 then 
    Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
  fi;
  if isgeom=true then
    Printn("kk", kk, "isgeom");
    if NoOverlattice(tGram, th)<>true then buzz(476127461); fi;
    scons:=GetTotalSplcons(tGram, th);
    if Trec.count=0 then buzz(47652);fi;
    if Trec.count>1 then Printn("____ more than one T", Trec.count); fi;
    trec.scons:=scons;
    trec.Trec:=Trec;
    Printn("____ scons", Length(scons)); 
  else
    if tsign[3]>0 and Trec.count<>0 then buzz(222652);fi;
    Printn("kk", kk, "not isgeom", isgeom);
  fi;
  #
  Add(AllZeroAdjRecs, trec);
od;

savedata(AllZeroAdjRecs);


# gap> Read("Task20260923Ue.g");
# All-zero adj: 
# kk 2 isgeom 
# ____ scons 3 
# kk 3 isgeom 
# ____ scons 6 
# kk 4 isgeom 
# ____ scons 10 
# kk 5 isgeom 
# ____ scons 15 
# kk 6 isgeom 
# ____ scons 21 
# kk 7 isgeom 
# ____ scons 28 
# kk 8 isgeom 
# ____ scons 36 
# kk 9 isgeom 
# ____ scons 45 
# kk 10 isgeom 
# ____ scons 55 
# kk 11 not isgeom [ false, 4 ] 
# kk 12 not isgeom [ false, 4 ] 
# kk 13 not isgeom [ false, 4 ] 
# kk 14 not isgeom [ false, 4 ] 
# kk 15 not isgeom [ false, 4 ] 
# kk 16 not isgeom [ false, 4 ] 
# kk 17 not isgeom [ false, 4 ] 
# kk 18 not isgeom [ false, 4 ] 
# kk 19 not isgeom [ false, 4 ] 






############