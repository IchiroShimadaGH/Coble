#Read("CheckPaper11.g");

Read("paperaffineConway/affineConwayCompdata.txt");
Read("AutDiscfByStabilizerChain.txt");
Read("checktools.g");





tG:=erecs[3].minusrec.Gram;
ttG:=2*InverseMat(tG);
SignatureQ(ttG);
if ShortestVectors(ttG, 1).vectors<>[] then beep(81881); fi;


Printn("Check the norbit decomp by stabilizer");
Printn("_______________ 2026/08/20/12:00");
