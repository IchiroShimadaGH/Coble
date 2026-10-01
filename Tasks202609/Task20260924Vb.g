#Read("Task20260924Vb.g");

Read("SplitConsTools.g");
Read("even_lattice_genus.g");
Read("NonSingOverLats.g");



#bb:=2;

totalrecsname:=Concatenation("totalrecsname", String(bb));
Printn(totalrecsname);

TotalRecs:=[];

for kk in [2..19] do 
  Printn("_______________");
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
  oreccounter:=0;
  for orec in orecs do
    oreccounter:=oreccounter+1;
    ttrec:=StructuralCopy(trec);
    ttrec.orec:=orec;
    tGram:=orec.Gram;
    th:=orec.h;
    isgeom:=IsGeom(tGram, th);
    ttrec.isgemo:=isgeom;
    # ttmdiscrec:=DiscriminantForm(-tGram);
    # tsign:=[22, 3, 19]-SignatureQ(tGram);
    # ttsign:=[tsign[2], tsign[3]];
    # if tsign[3]>0 then 
    #   Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
    # fi;
    if isgeom=true then
      sprats:=GetSpRats(tGram, th);
      # if Trec.count=0 then buzz(47652);fi;
      # if Trec.count>1 then Printn("____ more than one T", Trec.count); fi;
      ttrec.sprats:=sprats;
      # trec.Trec:=Trec;
      Printn("bb", bb, "kk", kk, "isgeom", 
          "extdeg", orec.extdeg, "sprats", sprats.nopss,   
          ":",  oreccounter, "in", Length(orecs)); 
    else
      # if tsign[3]>0 and .count<>0 then buzz(222652);fi;
      Printn("bb", bb, "kk", kk, "not isgeom", isgeom,
      ":",  oreccounter, "in", Length(orecs));
    fi;
    #
    Add(TotalRecs, ttrec);
  od;
od;

savedataas(TotalRecs, totalrecsname);




#######################