
#Read("Task20261002Ub.g");

readdata("FoundShRecsOverLats");

Read("NonSingOverLats.g");
Read("even_lattice_genus.g");
Read("SplitConsTools.g");
Read("ReidemeisterSchreier.g");
Read("AutShByGramP.g");

# gap> RecNames(FoundShRecsOverLats[1]);
# [ "h", "rho", "Gram", "basisdual", "extgroupgens", "extdeg", "ampleflag", "Trec", "sprats" ]

theShrecs:=[];
tname:="temp2Shrecs";



counter:=0;
for trec in FoundShRecsOverLats do
  counter:=counter+1;
  if trec.rho<>Length(trec.Gram) then Error("Inconsistent rho"); fi;
  if trec.rho>7 then continue; fi;
  if trec.Trec.count<>1 then Error("Expected a unique T"); fi;
  if trec.sprats.scons=[] or RankMat(Concatenation(trec.sprats.scons))<trec.rho then 
    continue; 
  fi;
  IsNewPic(trec);
  Printn("____", counter, "in" , Length(FoundShRecsOverLats));
od;

savedataas(theShrecs, tname);




##########

