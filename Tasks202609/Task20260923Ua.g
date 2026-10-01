

#Read("Task20260923Ua.g");

Read("SplitConsTools.g");
Read("even_lattice_genus.g");

readdata("randomoddenhancedCopy");
readdata("randomoddenhanced2Copy");

max1:=Maximum(List(randomoddenhanced, xx->xx[1]));
Printn(max1);

poss:=Positions(List(randomoddenhanced, xx->xx[1]), max1);
maxadj1:=randomoddenhanced[poss[1]][2];;

trec:=WGraphToGramh(maxadj1, Length(maxadj1));;
tGram:=trec.Gram;
th:=trec.h;
isgeom:=IsGeom(tGram, th);
if isgeom<>true then beep(27162); fi;

scons:=GetTotalSplcons(tGram, th);
if Length(scons)<>max1 then beep(9192992); fi;


ttmdiscrec:=DiscriminantForm(-tGram);
tsign:=[22, 3, 19]-SignatureQ(tGram);
ttsign:=[tsign[2], tsign[3]];

 Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
if Trec.count=0 then buzz(47652);fi;

trec.scons:=scons;
trec.Trec:=Trec;

savedataas(trec, "Data1845");


############