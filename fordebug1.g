#readdata("sublat180rec");

tGramS:=DiagonalMats([[[8]], -sublat180rec.Gram]);
tdiscrec:=DiscriminantForm(tGramS);
periodcond:=[IdentityMat(2), DiagonalMat(List(tdiscrec.discg, xx->xx-1))];

geomauts:=[];


for tg in sublat180rec.allgs do 
  ttg:=DiagonalMats([[[1]], tg]);
  tqg:=OLtoOqL(ttg, tdiscrec);
  if tqg in periodcond then 
    Add(geomauts, tg);
  fi;
od;

Printn(Length(geomauts));

