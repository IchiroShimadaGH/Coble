#Read("AutShByGramP.g");


GetGramP:=function(GramS, h)
  local hdual, GramP;
  hdual:=h*GramS;
  GramP:=TransposedMat([hdual])*[hdual]-GramS;
  if h*GramP*h<>2 then beep(998912); fi;
  return(GramP);
end;

AutShByGramP:=function(GramS, h)
  local n, GramP, basisrecP, ogrec, gs, hs, pos, th, tgen, thtgen, 
  tgenspermsonhs, tg, tgp, hstabwds,  hstabgens, twd, xx,  ttg, tth;
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
  for twd in hstabwds do 
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
    AutShsize:=ogrec.size/Length(hs)
  );
  return(trec);
end;

IsIsomShsByGramP:=function(Sh1rec, GramS2, h2, basisrecP2)
end;




#############