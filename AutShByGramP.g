#Read("AutShByGramP.g");

#Read("ReidemeisterSchreier.g");


DiscfImageIndex:=function(OGgens, GramL) 
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
  if leng=0 then 
    return rec(sizeaut:=1,sizeimage:=1,indeximage:=1); 
  fi;
  vs:=Set(Cartesian(List(discg, aa->[0..aa-1])));
  bigmat:=CopyNormalDiscf(TMTTmult(vs, discrec.discf));
  inibiis:=List(IdentityMat(leng), tv->SinglePosition(vs, tv));
  OqLrec:=AutDiscfByStabilizerChain(inibiis, bigmat);
  sizeaut:=OqLrec.order;
  tpermgens:=[];
  for tg in OGgens do 
    tqg:=OLtoOqL(tg, discrec);
    vstqg:=List(vs*tqg, ttv->List([1..leng], ii-> (ttv[ii] mod discg[ii])));
    poss:=List(vstqg, ttv->PositionSorted(vs, ttv));
    Add(tpermgens, PermList(poss));
  od;
  sizeimage:=Size(Group(tpermgens,()));
  indeximage:=sizeaut/sizeimage;
  if not IsInt(indeximage) then beep(615825885); fi;
  issuujrec:=rec(
    sizeaut:=OqLrec.order, 
    sizeimage:=sizeimage,
    indeximage:=indeximage
  );
  return(issuujrec);
end;


GetGramP:=function(GramS, h)
  local hdual, GramP;
  hdual:=h*GramS;
  GramP:=TransposedMat([hdual])*[hdual]-GramS;
  if h*GramP*h<>2 then beep(998912); fi;
  return(GramP);
end;

AutShByGramP:=function(GramS, h)
  local n, GramP, basisrecP, ogrec, gs, hs, pos, th, tgen, thtgen, 
  tgenspermsonhs, tg, tgp, hstabwds,  hstabgens, twd, xx,  ttg, tth,
  etarec, trec;
  #
  n:=Length(h);
  GramP:=GetGramP(GramS, h);
  basisrecP:=BasisRec(GramP);
  ogrec:=OGLatFromBasisRec(basisrecP);
  gs:=[IdentityMat(n)];
  hs:=[List(h)];
  pos:=0;
  for th in hs do
    pos:=pos+1;
    for tgen in ogrec.totalgens do 
      thtgen:=th*tgen;
      if not thtgen in hs then 
        Add(hs, thtgen);
        Add(gs, gs[pos]*tgen);
      fi;
    od;
  od;
  #
  tgenspermsonhs:=[];
  for tg in ogrec.totalgens do 
    tgp:=PermList(List(hs*tg, tv->SinglePosition(hs, tv)));
    Add(tgenspermsonhs, tgp);
  od;
  #
  hstabwds:=ReidemeisterSchreierStabilizer(tgenspermsonhs);
  hstabgens:=[];
  #
  for twd in hstabwds.integerWords do 
    ttg:=IdentityMat(n);
    for xx in twd do 
      if xx>0 then ttg:=ttg*ogrec.totalgens[xx];
      elif xx<0 then 
        tth:=InverseMat(ogrec.totalgens[-xx]);
        ttg:=ttg*tth;
      else buzz(576576); 
      fi;
    od;
    if h*ttg<>h then buzz(5157641); fi;
    Add(hstabgens, ttg);
  od;
  #
  #etarec:=DiscfImageIndex(hstabgens, GramS);
  #
 trec:=rec(
    rank:=Length(h),
    GramS:=GramS,
    h:=h,
    GramP:=GramP,
    basisrecP:=basisrecP, 
    ogrecP:=ogrec,
    horbit:=hs,
    transporters:=gs,
    AutShgens:= hstabgens,
    AutShsize:=ogrec.size/Length(hs),
    #etarec:=etarec
  );
  return(trec);
end;

IsIsomShsByGramP:=function(Sh1rec, GramS2, h2, basisrecP2)
  local beep, tg, tginv, h2tginv, pos, Ag, ttg;
  #
  beep:=function(beepnumb)
    localbeep("IsIsomShsByGramP", beepnumb); Error();
  end;
  #
  #
  #IsIsomBasisRecs:=function(basisrecA, basisrecB)
  tg:=IsIsomBasisRecs(Sh1rec.basisrecP, basisrecP2);
  if tg=false then return(false); fi;
  if TMTTmult(tg, basisrecP2.Gram)<>Sh1rec.GramP then beep(423142); fi;
  tginv:=InverseMat(tg);
  h2tginv:=h2*tginv;
  if h2tginv*Sh1rec.GramP*h2tginv<>2 then beep(561321); fi;
  if not h2tginv in Sh1rec.horbit then return(false); fi;
  pos:=SinglePosition(Sh1rec.horbit, h2tginv);
  Ag:=Sh1rec.transporters[pos];
  ttg:=Ag*tg;
  if TMTTmult(ttg, GramS2)<>Sh1rec.GramS then beep(49942); fi;
  if Sh1rec.h*ttg<>h2 then beep(524131); fi;
  return(ttg);
end;


theShrecs:=[];
tname:="temp2Shrecs";


IsNewPic:=function(olrec)
  #
  local GramS2, h2, basisrecP2, isnewflag, oldShrec, 
  ttg, newrec, discgT;
  #
  GramS2:=olrec.Gram;
  h2:=olrec.h;
  basisrecP2:=BasisRec(GetGramP(GramS2, h2));
  isnewflag:=true;
  for oldShrec in theShrecs do 
    if oldShrec.rank<>Length(GramS2) then continue; fi;
    if oldShrec.detS<>AbsInt(DeterminantIntMat(GramS2)) then continue;  fi;
    if oldShrec.sprats.nopss<>olrec.sprats.nopss then continue;  fi;
    ttg:=IsIsomShsByGramP(oldShrec, GramS2, h2, basisrecP2);
    if ttg<>false then 
      isnewflag:=false;
      Printn("__________________________________isom");
      break; 
    fi;
  od;
  if isnewflag then 
    newrec:=AutShByGramP(GramS2, h2);
    newrec.detS:=AbsInt(DeterminantIntMat(GramS2));
    newrec.sprats:=olrec.sprats;
    #newrec.Trec:=olrec.Trec;
    discgT:=DiscriminantForm(-GramS2).discg;
    if 22-Length(GramS2)<Length(discgT)+2 then 
      buzz(5616751); 
    else 
      newrec.Nikulin:=true;
    fi;
    Add(theShrecs, newrec);
    Printn("newrec", Length(theShrecs));
    savedataas(theShrecs, tname);
    # if newrec.etarec.indeximage>1 then  
    #   Add(thespecialShrecs, newrec);
    #   Printn("special newrec", Length(thespecialShrecs));
    #   savedataas(thespecialShrecs, t2name);
    # fi;
  fi;
  return();
end;




#############