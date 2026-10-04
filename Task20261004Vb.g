#Read("Task20261004Vb.g");


Read("Nipps.g");

NtAs:=[];


for tA in NippGramMatrices do
  diags:=List([1..5], ii->tA[ii][ii]);
  #if Set(diags)<>[4]  then continue; fi;
  if 2 in diags then continue; fi;
  detA:=DeterminantIntMat(tA);
  if detA>500 then continue; fi;
  sntA:=SmithNormalFormIntegerMat(tA);
  diagssntA:=List([1..5], ii->sntA[ii][ii]);
  notones:=Filtered(diagssntA, xx-> xx<>1);
  if Length(notones)>=3 then continue; fi;
  evens:=Filtered(diagssntA, xx-> xx mod 2=0);
  if Length(evens)>1 then continue; fi;
  Printn(DiscriminantForm(tA).discg, detA);
  Add(NtAs, tA);
od;



