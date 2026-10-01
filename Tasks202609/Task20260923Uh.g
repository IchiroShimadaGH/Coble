#Read("Task20260923Uh.g");

Read("SplitConsTools.g");
Read("even_lattice_genus.g");
Read("NonSingOverLats.g");



bb:=0;

totalrecsname:=Concatenation("totalrecsname", String(bb));
Printn(totalrecsname);

TotalRecs:=[];

for kk in [2..19] do 
  adjmat:=(-2)*IdentityMat(kk);
  for ii in [1..kk] do 
    for jj in [ii+1..kk] do 
      adjmat[ii][jj]:=bb;
      adjmat[jj][ii]:=bb;
    od;
  od;
  trec:=WGraphToGramh(adjmat, kk);
  t0Gram:=trec.Gram;
  t0h:=trec.h;
  orecs:=NonSingOverLats(t0Gram,  t0h);
  Printn("kk", kk, "orecs", Length(orecs));
  for orec in orecs do
    tGram:=orec.Gram;
    th:=orec.h;
    isgeom:=IsGeom(tGram, th);
    trec.isgemo:=isgeom;
    # ttmdiscrec:=DiscriminantForm(-tGram);
    # tsign:=[22, 3, 19]-SignatureQ(tGram);
    # ttsign:=[tsign[2], tsign[3]];
    # if tsign[3]>0 then 
    #   Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
    # fi;
    if isgeom=true then
      Printn("____kk", kk, "isgeom");
      sprats:=GetSpRats(tGram, th);
      # if Trec.count=0 then buzz(47652);fi;
      # if Trec.count>1 then Printn("____ more than one T", Trec.count); fi;
      trec.sprats:=sprats;
      # trec.Trec:=Trec;
      Printn("____ sprats", sprats.nopss); 
    else
      # if tsign[3]>0 and .count<>0 then buzz(222652);fi;
      Printn("____kk", kk, "not isgeom", isgeom);
    fi;
    #
    Add(TotalRecs, trec);
  od;
od;

savedataas(TotalRecs, totalrecsname);




#######################