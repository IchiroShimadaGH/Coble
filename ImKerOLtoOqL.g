#Read("ImKerOLtoOqL.g");




ImOLtoOqL:=function(OGgens, GramL)
  local discrec, discg, vs, aa, leng, tpermgens, 
  tg, tqg, vstqg, ii, poss, size, ttv;
  discrec:=DiscriminantForm(GramL);
  discg:=discrec.discg;
  vs:=Set(Cartesian(List(discg, aa->[0..aa-1])));
  leng:=Length(discg);
  tpermgens:=[];
  for tg in OGgens do 
    tqg:=OLtoOqL(tg, discrec);
    vstqg:=List(vs*tqg, ttv->List([1..leng], ii-> (ttv[ii] mod discg[ii])));
    poss:=List(vstqg, ttv->Position(vs, ttv));
    Add(tpermgens, PermList(poss));
  od;
  size:=Size(Group(tpermgens));
  return(size);
end;

KerOLtoOqL:=function(GramL)
  local beep, nn, GramLdual, dd, ddGramLdual, LLLrec, thebs, bsnorms, 
  svrec, nrmsset, posss, pvss, getpvs, poss, pos, isinL, 
  candidatess, candidatesduals, tb, candidates, pvs, tpv, tv, candidatesdual, 
  newGram, intnumbss, jj, ii, thebsinv, discrec, idqg, reachflag, thetg, 
  simplexxt, leng, tgdual, tg, tcans, tcansdual, intnumbs, tvpos, tu, stabrecs,
   tsize, tgens, tlevel, gens, nops, initbs, tbs, initbsdual, iniparsol, 
   stabrec, kerrec;
  #
  beep:=function(beepnumb)
    localbeep("KerOLtoOqL", beepnumb); Error();
  end;
  #
  nn:=Length(GramL);
  GramLdual:=InverseMat(GramL);
  dd:=Lcm(List(Flat( GramLdual), DenominatorRat));
  ddGramLdual:=dd*GramLdual;
  LLLrec:=LLLReducedGramMat(ddGramLdual);
  thebs:=LLLrec.transformation;
  bsnorms:=List(thebs, tb->tb*ddGramLdual*tb);
  svrec:=ShortestVectors(ddGramLdual, Maximum(bsnorms));
  nrmsset:=Set(svrec.norms);
  posss:=List(nrmsset, aa->Positions(svrec.norms, aa));
  pvss:=List(posss, poss->List(poss, pos->svrec.vectors[pos]));
  #
  getpvs:=function(aa)
    return(pvss[SinglePosition(nrmsset, aa)]);
  end;
  #
  isinL:=function(vv)
    return(IsZeroVect((vv*ddGramLdual) mod dd));
  end;
  #
  candidatess:=[];
  candidatesduals:=[];
  for tb in thebs do
    candidates:=[];
    pvs:=getpvs(tb*ddGramLdual*tb);
    for tpv in pvs do 
      for tv in [tpv, -tpv] do 
        if isinL(tv-tb) then Add(candidates, tv); fi;
      od;
    od;
    candidatesdual:=candidates*ddGramLdual;
    Add(candidatess, candidates);
    Add(candidatesduals, candidatesdual); 
  od;
  #
  newGram:=TMTTmult(thebs, ddGramLdual);
  intnumbss:=List([1..nn], jj->List([1..jj-1], ii->newGram[jj][ii]));
  thebsinv:=InverseMat(thebs);
  discrec:=DiscriminantForm(GramL);
  idqg:=IdentityMat(Length(discrec.discg));
  reachflag:=false;
  thetg:=[];
  #
  simplexxt:=function(parsol)
    local tleng, tgdual, tg, tu, tcans, tcansdual, 
    intnumbs, tvpos, tv;
    tleng:=Length(parsol);
    if tleng=nn then 
      reachflag:=true;
      tgdual:=thebsinv*parsol;
      tg:=GramL*tgdual*GramLdual;
      if not IsIntMat(tg) then beep(58811); fi;
      if TMTTmult(tg, GramL)<>GramL then beep(881811); fi;
      if OLtoOqL(tg, discrec)<>idqg then beep(9192362); fi;
      thetg:=ShallowCopy(tg);
      return();
    fi;
    tcans:=candidatess[tleng+1];
    tcansdual:=candidatesduals[tleng+1];
    intnumbs:=intnumbss[tleng+1];
    tvpos:=0;
    for tv in tcansdual do 
      tvpos:=tvpos+1;
      if parsol*tv =intnumbs  then 
        tu:=tcans[tvpos];
        Add(parsol, tu);
        simplexxt(parsol);
        if tu<>Remove(parsol) then beep(8127172); fi;
        if reachflag then return(); fi;
      fi;
    od;
    return();
  end;
  #
  stabrecs:=[];
  tsize:=1;
  tgens:=[];
  #
  #
  for tlevel in [1..nn] do 
    #Printn("start tlevel", tlevel, nn);
    gens:=[];
    nops:=0;
    tcans:=candidatess[tlevel];
    initbs:=List([1..tlevel-1], ii->thebs[ii]);
    if tlevel>1 then 
      initbsdual:=initbs*ddGramLdual;
      intnumbs:=intnumbss[tlevel];
    fi;
    for tv in tcans do 
      if tlevel=1 or initbsdual*tv=intnumbs then 
        iniparsol:=CopyAdd(initbs, tv);
        reachflag:=false;
        simplexxt(iniparsol);
        if reachflag then
          Add(gens, thetg);
          nops:=nops+1; 
        fi;
      fi;
    od;
    stabrec:=rec(level:=tlevel, gens:=gens, nops:=nops);
    Add(stabrecs, stabrec);
    tsize:=tsize*nops;
    tgens:=Union(tgens, gens);
  od;
  #
  kerrec:=rec(
    gens:=tgens, 
    size:=tsize
  );
  return(kerrec);
end;

CheckKerq:=function(GramL, gs)
  local beep, svrec, pvs, vs, permgens, tnrm,
  tdiscrec, idleng, tpg, tg, tbasis, tv;
  #
  tdiscrec:=DiscriminantForm(GramL);
  idleng:=IdentityMat(Length(tdiscrec.discg));
  for tg in gs do 
    if OLtoOqL(tg, tdiscrec)<>idleng then beep(583211); fi;
  od;
  return(true);
  #
end;

################
