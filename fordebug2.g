
 GramLをGram行列にもつ偶格子 Lの自己同型群の部分集合OGgensがあたえられたとき，
 OGgensで生成される O(L) の部分群 G の自然な準同型
 O(L) to O(qL) (qLは Lのdiscriminant form ね)
 による像のindex を求める関数 IsSurjOLtoOqL を書いたの.
 sizeautはO(qL) のサイズ，
 sizeimage Gの像のサイズね．
 批評して修正点を指摘して：

IsSurjOLtoOqL:=function(OGgens, GramL) 
  local discrec, discg, vs, aa, leng, tpermgens, 
  tg, tqg, vstqg, ii, poss, sizeaut, ttv, sizeimage,
  indeximage, bigmat, inibiis, OqLrec, issuujrec, beep;
  #
  beep:=function(beepnumb)
    localbeep("IsSurjOLtoOqL", beepnumb); Error();
  end;
  #
  discrec:=DiscriminantForm(GramL);
  discg:=discrec.discg;
  leng:=Length(discg);
  vs:=Set(Cartesian(List(discg, aa->[0..aa-1])));
  bigmat:=CopyNormalDiscf(TMTTmult(vs, discrec.discf));
  inibiis:=List(IdentityMat(leng), tv->SinglePosition(vs, tv));
  OqLrec:=AutDiscfByStabilizerChain(inibiis, bigmat);
  sizeaut:=OqLrec.order;
  #
  tpermgens:=[];
  for tg in OGgens do 
    tqg:=OLtoOqL(tg, discrec);
    vstqg:=List(vs*tqg, ttv->List([1..leng], ii-> (ttv[ii] mod discg[ii])));
    poss:=List(vstqg, ttv->Position(vs, ttv));
    Add(tpermgens, PermList(poss));
  od;
  sizeimage:=Size(Group(tpermgens));
  indeximage:=sizeaut/sizeimage;
  if not IsInt(indeximage) then beep(615825885); fi;
  issuujrec:=rec(
    sizeaut:=sizeaut, 
    sizeimage:=sizeimage,
    indeximage:=indeximage
  );
  return(issuujrec);
end;