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


#theaut:=SimpleOGLat(orthorec.Gram);;
# これは時間がかかりすぎ．
# ChatGPTを経由して Pariで計算してもらった
#

# 2 0 [ [ [ 12, 0 ], [ 0, 30 ] ] ] 
# 4 0 [  ] 
# 8 0 [ [ [ 12, 0 ], [ 0, 120 ] ] ] 
# 10 0 [  ] 
# 14 0 [ [ [ 30, 0 ], [ 0, 84 ] ] ] 
# 16 0 [  ] 
# 20 0 [ [ [ 60, 0 ], [ 0, 60 ] ] ] 
# 22 0 [  ] 
# 26 0 [ [ [ 6, 0 ], [ 0, 780 ] ] ] 
# 28 0 [  ] 
# 32 0 [ [ [ 12, 0 ], [ 0, 480 ] ] ] 
# 34 0 [  ] 
# 38 0 [ [ [ 12, 0 ], [ 0, 570 ] ] ] 
# 40 0 [  ] 
# 44 0 [ [ [ 84, 12 ], [ 12, 96 ] ] ] 
# 46 0 [  ] 
# 50 0 [ [ [ 30, 0 ], [ 0, 300 ] ] ] 
# 52 0 [  ] 
# 56 0 [ [ [ 84, 0 ], [ 0, 120 ] ] ] 
# 58 0 [  ] 
# 62 0 [ [ [ 12, 0 ], [ 0, 930 ] ], [ [ 30, 0 ], [ 0, 372 ] ] ] 
# 64 0 [  ] 
# 68 0 [ [ [ 12, 0 ], [ 0, 1020 ] ] ] 


# sublat 180から得られるdegree 4, 6のK3は存在しない．
