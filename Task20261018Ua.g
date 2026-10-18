
 #Read("Task20261018Ua.g");

#readdata("new2singK3configs");

#  gap>  Set(new2singK3configs, RecNames);
# [ [ "h", "GramS", "rr", "sprats", "tt", "dd", 
# "GammaSize", "LambdaSize", "gtypes", "hperprec", "nonsingOLs" ] ]

#  Set(new2singK3configs, xx->RecNames(xx.hperprec));

# gap> Set(new2singK3configs, xx->RecNames(xx.hperprec));
# [ [ "h", "T", "rank", "GramS", "AutShgens", "AutShsize", "GramN", "Tinv", "Nbasisrec" ] ]

the4455rec:=new2singK3configs[1];;

discSXrec:=DiscriminantForm(the4455rec.GramS);
discTXrec:=DiscriminantForm((-1)*the4455rec.GramS);

TXcandidates:=EvenPosDefBinaryFormGL(Product(discTXrec.discg));; 

#GeneralizedAnIsomFQF:=function(discf1, discf2)

the4455rec.discmSXrec:=DiscriminantForm((-1)*the4455rec.GramS);

discfmSX:=[the4455rec.discmSXrec.discg, the4455rec.discmSXrec.discf];

for tTX in TXcandidates do 
  disctTXrec:=DiscriminantForm(tTX);
  discfTX:=[disctTXrec.discg, disctTXrec.discf];
  anisom:=GeneralizedAnIsomFQF(discfmSX,  discfTX);
  if anisom<>fail then 
    if CopyNormalDiscf(TMTTmult(anisom, disctTXrec.discf))<>discfmSX[2] then beep(696111); fi;
    Printn("found", tTX); 
    the4455rec.GramTX:=tTX;
    the4455rec.discTXrec:=disctTXrec;
    the4455rec.anisomSXTX:=List(anisom);
  fi;
od;


###############