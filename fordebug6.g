#Read("Task20261025Ua.g");



#readdata("singK3configs");

tt:=singK3configs[1].tt;

ttvs:=List(tt, ii->MakeVectei(24, ii));

#readdata("LeechBasis");
readdata("GramLeech");

orthorec:=OrthogonalCompRec(GramLeech, ttvs);

Determinant(orthorec.Gram);

#GeneralizedAnIsomFQF:=function(discf1, discf2)

for td in 2*[1..200] do 
  tGram:=DiagonalMats([ [[td]] , -orthorec.Gram] );;
  tdiscrec:=DiscriminantForm(tGram);
  if Length(tdiscrec.discg)=2 then 
    tdiscf:=[tdiscrec.discg, tdiscrec.discf];
    tdet:=Product(tdiscrec.discg);
    if not tdet mod 4 in [3, 0] then continue; fi;
    tTs:=EvenPosDefBinaryFormGL(tdet);
    geomTs:=[];
    for GramtT in EvenPosDefBinaryFormGL(tdet) do 
      discrecmT:=DiscriminantForm(-GramtT);
      tdiscf2:=[discrecmT.discg, discrecmT.discf];
      #GeneralizedAnIsomFQF:=function(discf1, discf2)
      anisom:=GeneralizedAnIsomFQF( tdiscf,  tdiscf2);
      if anisom<>fail then 
        Add(geomTs, GramtT);
      fi;
    od;
    Printn(td,  geomTs);
  else 
    Printn(td, tdiscrec.discg);
  fi;
od;


theaut:=SimpleOGLat(orthorec.Gram);;