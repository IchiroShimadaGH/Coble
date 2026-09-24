#Read("Task20260924Ue.g");

Read("SplitConsTools.g");
Read("even_lattice_genus.g");
Read("NonSingOverLats.g");

readdata("Data1845");

adj1845:=Data1845.adj;

kk:=Length(adj1845);  #=17

TT:=[[-2]];

OverOL1845:=[];

counter:=0;

iter:=IteratorOfCartesianProduct(List([1..17], ii->[1,3]));
for tAv in iter do 
  counter:=counter+1;
  if counter mod 1000 =0 then Printn("___", counter, 2^17); fi;
  tA:=[tAv];
  ttA:=TransposedMat(tA);
  newadj:=MatMatToMat([[adj1845, ttA], [tA, TT]]);
  trec:=WGraphToGramh(newadj, kk+1);
  tGram:=trec.Gram;
  th:=trec.h;
  tsign:=SignatureQ(tGram);
  if tsign[2]<>1 then continue; fi;
  orecs:=NonSingOverLats(tGram, th); 
  Printn("signature", tsign, "orecs", Length(orecs));
  for orec in orecs do 
    ttGram:=orec.Gram;
    tth:=orec.h;
    isgeom:=IsGeom(ttGram, tth);
    if isgeom=true then 
      sprats:=GetSpRats(ttGram, tth);
      ttrec:=rec(
        adj:=newadj, 
        iniGram:=tGram,
        Gram:=ttGram,
        inih:=th, 
        h:=tth, 
        sprats:=sprats,
        addwds:=orec.addwds
      );
      Add(OverOL1845, ttrec);
      savedata(OverOL1845);
      Printn(orec.extdeg, sprats.nopss, Length(OverOL1845));
    else 
      Printn(isgeom, orec.extdeg, sprats.nopss);
    fi;
  od;#for orec in orecs do 
od;









############