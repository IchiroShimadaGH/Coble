
#Read("Task20260923Ub.g");

Read("SplitConsTools.g");
Read("even_lattice_genus.g");

readdata("randomoddenhancedCopy");
readdata("randomoddenhanced2Copy");

randomodds:=Concatenation(randomoddenhanced, randomoddenhanced2);
totalleng:=Length(randomodds);
Printn(totalleng);

RandomOdds:=[];

counter:=0;
for tdata in randomodds do 
  counter:=counter+1;
  nopsscons:=tdata[1];
  tadj:=tdata[2];
  trec:=WGraphToGramh(tadj, Length(tadj));;
  tGram:=trec.Gram;
  th:=trec.h;
  isgeom:=IsGeom(tGram, th);
  if isgeom<>true then beep(27162); fi;
  scons:=GetTotalSplcons(tGram, th);
  if Length(scons)<>nopsscons then beep(9192992); fi;
  ttmdiscrec:=DiscriminantForm(-tGram);
  tsign:=[22, 3, 19]-SignatureQ(tGram);
  ttsign:=[tsign[2], tsign[3]];
  Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
  if Trec.count=0 then buzz(47652);fi;
  trec.scons:=scons;
  trec.Trec:=Trec;
  Add(RandomOdds, trec);
  savedataas(RandomOdds, "RandomOdds20260923");
  Printn(counter,"in", totalleng, ":", Length(tGram), nopsscons, ttmdiscrec.discg);
od;








############