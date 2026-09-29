#Read("OGLat.g");

#
# made by codex on 2026/09/28
#

VssRec:=function(GramL, nrmcandidate)
  #
  # Enumerate all nonzero vectors up to sign with norm <= nrmcandidate.
  # Require that they generate L over Z. minreqnrm is the first
  # generating shell; minflag says that the bound equals that norm.
  #
  local nn, svsrec, nrmsset, vects, vss, tvs, tnrm, 
  tpos, doesgenerate, vssrec, poss, vs, minflag,
  minreqnrm; #minimal required norm 
  #
  nn:=Length(GramL);
  svsrec:=ShortestVectors(GramL, nrmcandidate);
  nrmsset:=Set(svsrec.norms);
  if nrmsset=[] then Error("VssRec: no nonzero vectors within bound"); fi;
  vects:=svsrec.vectors;
  vss:=[];
  tvs:=[];
  minflag:=true; 
  minreqnrm:=fail;
  for tnrm  in nrmsset do 
    poss:=Positions(svsrec.norms, tnrm);
    vs:=List(poss, tpos->vects[tpos]);
    Add(vss, vs);
    Append(tvs, vs);
    doesgenerate:=(Rank(tvs)=nn and CokerTorsion(tvs)=[]);
    if doesgenerate and minreqnrm=fail then
      minreqnrm:=tnrm;
    fi;
  od;
  if minreqnrm=fail then
    Error("VssRec: short vectors do not generate the lattice over Z");
  fi;
  minflag:=(minreqnrm=nrmcandidate);
  vssrec:=rec(
    Gram:=List(GramL), 
    maxnrm:=nrmcandidate,
    nrmsset:=nrmsset, 
    nopss:=List(vss, Length),
    minflag:=minflag, 
    minreqnrm:=minreqnrm,
    vss:=vss
  );
  return(vssrec);
end;


Webrint:=function(tv, Gram, tvs)
  local tvdual;
  tvdual:=tv*Gram;
  return(Collected(List(tvs*tvdual, AbsInt)));
end;

Webrints:=function(tv, Gram, tvss)
  local tvs;
  return(List(tvss, tvs->Webrint(tv, Gram, tvs)));
end;


BasisRec:=function(arg)
  #
  local beep, GramL, trialnumb, nn, norm, tU, tt,GramL2,  LLLrec, basis, ii,  
  basisnrms,fflag,  maxnrm,  tbasisnrms, 
  newGram, getvs, trec, jj, kk, tintnumbs,
  tbasisdual, tv, tvs,  tbasis, thebasis,
  nrmcandidate, ttintnumbs, maxbasisrrms,  tmaxbasisnrms,
  counter, minreqnrm, thevss,  basiswebrintss,
  candidatess, bpos, tcandidates, inipos, subcandidatess, tsubcandidates,
  tbintnumbs, fingerprints, candidatessleng, perm;
  #
  #
  beep:=function(beepnumb)
    localbeep("BasisRec", beepnumb); Error();
  end;
  #
  if Length(arg)=1 then 
    GramL:=arg[1];
    trialnumb:=10;
  elif Length(arg)=2 then 
    GramL:=arg[1];
    trialnumb:=arg[2];
  else beep(3112221);
  fi;
  #
  nn:=Length(GramL);
  if SignatureQ(GramL)<>[nn, nn, 0] then beep(665522); fi;
  fflag:=false;
  basisnrms:=[];
  maxbasisrrms:=infinity;
  #
  for tt in [1..5] do
    tU:=RandomUnimodMat(nn);
    GramL2:=TMTTmult(tU, GramL);
    LLLrec:=LLLReducedGramMat(GramL2);
    tbasis:=LLLrec.transformation;
    tbasisnrms:=List([1..nn], ii->LLLrec.remainder[ii][ii]);
    SortParallel(tbasisnrms, tbasis);
    tmaxbasisnrms:=Maximum(tbasisnrms);
    if basisnrms=[] or 
      maxbasisrrms> tmaxbasisnrms or 
      (maxbasisrrms= tmaxbasisnrms and basisnrms>tbasisnrms) then
      thebasis:=tbasis*tU;
      basisnrms:=List(tbasisnrms);
      maxbasisrrms:= tmaxbasisnrms;
    fi;
  od;
  #
  nrmcandidate:=maxbasisrrms;
  trec:=VssRec(GramL, nrmcandidate);
  #
  if trec.minflag=false then 
    minreqnrm:=trec.minreqnrm;
    fflag:=false;
    counter:=0;
    while not fflag and counter<trialnumb do
      counter:=counter+1;
      tU:=RandomUnimodMat(nn);
      GramL2:=TMTTmult(tU, GramL);
      LLLrec:=LLLReducedGramMat(GramL2);
      tbasisnrms:=List([1..nn], ii->LLLrec.remainder[ii][ii]);
      if Maximum(tbasisnrms)=minreqnrm then 
        fflag:=true;
        tbasis:=LLLrec.transformation;
        SortParallel(tbasisnrms, tbasis);
        thebasis:=tbasis*tU;
        basisnrms:=List(tbasisnrms);
        maxbasisrrms:= minreqnrm;
        trec:=VssRec(GramL, minreqnrm);
      fi;
    od;
  fi;
  #
  #
  thevss:=trec.vss;
  #
  #
  basiswebrintss:=List(thebasis, tb->Webrints(tb, GramL, thevss));
  #
  candidatess:=[];
  for bpos in [1..nn] do 
    tvs:=thevss[SinglePosition(trec.nrmsset, basisnrms[bpos])];
    tcandidates:=[];
    for tv in tvs do 
      if Webrints(tv, GramL, thevss)=basiswebrintss[bpos] then 
        Add(tcandidates, tv);
      fi;
    od;
    #
    Add(candidatess, tcandidates);
  od;
  #
  #
  candidatessleng:=List(candidatess, Length);
  perm := Sortex(candidatessleng);
  #
  candidatess:=Permuted(candidatess, perm);
  inipos:=function(tv)
    local xx;
    for xx in tv do 
      if xx>0 then return(tv); elif xx<0 then return(-tv); fi;
    od;
    beep(727211);
  end;
  candidatess[1]:=List(candidatess[1], inipos);
  # to make the computation in getpermcan1 easy
  trec.candidatess:=candidatess;
  #
  basisnrms:=Permuted(basisnrms, perm);
  trec.basisnrms:=basisnrms;
  thebasis:=Permuted(thebasis, perm);
  trec.basis:=thebasis;
  basiswebrintss:=Permuted(basiswebrintss, perm);
  trec.basiswebrintss:=basiswebrintss;
  newGram:=TMTTmult(thebasis, GramL);
  trec.newGram:=newGram;
  #
  #
  #
  ################
  #
  # A candidates is given by webrints 
  # A subcandidates is a subset defined by 
  # intersection numbers with younger vectors in the basis.
  #
  subcandidatess:=[];
  Add(subcandidatess, candidatess[1]);
  #
  for jj in [2..nn] do 
    tbintnumbs:=List([1..jj-1], kk->newGram[jj][kk]);
    tvs:=candidatess[jj];
    tbasisdual:=List([1..jj-1], kk->thebasis[kk])*GramL;
    tsubcandidates:=[];
    for tv in tvs do 
      ttintnumbs:=tbasisdual*tv;
      if ttintnumbs=tbintnumbs  then 
        Add(tsubcandidates, tv);
      fi;
      if ttintnumbs=-tbintnumbs  then 
        Add(tsubcandidates, -tv);
      fi;
      #
      # The above part cannot be changed to elif, 
      # because tintnumbs may be a zero vector
      #
    od;
    Add(subcandidatess, tsubcandidates);
  od;
  #
  fingerprints:=List(subcandidatess, Length);
  #
  trec.subcandidatess:=subcandidatess;
  trec.fingerprints:=fingerprints;
  #Printn(fingerprints);
  #
  return(trec);
end;


OGLatFromBasisRec:=function(arg) # arg is (basisrec) or (basisrec, giveupsec)
  #
 local beep, basisrec, giveupflag, giveupruntime, gvstopflag, st,
  GramL, nn, basis, newGram, candidatess, candidatessdual,
  subcandidatess, fingerprints, inipos, thevss, vpos, bpos,
  ttintnumbs, subcandidates,subcandidatessdual, 
  basisinv, basisdual, reachflag, thetg, simpleextend,
  subcandidates1, getpermcan1, stabrecs, doneposs,
  gengs, genperms, trueposs, neworb, tv, candidates,
  getpermcan, tperm, inipsol, tvduals, tvs, getorb,
  getorbsunion, pos1, falseseeds, trueseeds, falseposs,
  thelevel, inipsoldual, tbintnumbs, tpv, totalgens,
  lengcan, lengcan1, tlist, flist, pos, stabrec, OGrec;
  #
  beep:=function(beepnumb)
    localbeep("OGLatFromBasisRec", beepnumb); Error();
  end;
  #
  if Length(arg)=1 then 
    basisrec:=arg[1];
    giveupflag:=false;
    giveupruntime:=infinity;
  elif Length(arg)=2 then 
    basisrec:=arg[1];
    giveupflag:=true;
    giveupruntime:=1000*arg[2];
  else beep(212341);
  fi;
  #
  GramL:=basisrec.Gram;
  nn:=Length(GramL);
  basis:=basisrec.basis;
  newGram:=basisrec.newGram;
  if TMTTmult(basis, GramL)<>newGram then beep(391919); fi;
  basisinv:=InverseMat(basis);
  basisdual:=basis*GramL;
  #
  thevss:=basisrec.vss;
  #
  candidatess:=basisrec.candidatess;
  candidatessdual:=List(candidatess, tvs->tvs*GramL);
  subcandidatess:=basisrec.subcandidatess;
  subcandidatessdual:=List(subcandidatess, tvs->tvs*GramL);
  fingerprints:=basisrec.fingerprints;
  #
  gvstopflag:=false;
  reachflag:=false;
  thetg:=[];
  #
  simpleextend:=function(psol)
    local leng, tg, newsubcandidates,  tvs, 
    tpv, tv, tintnumbs, pos, ss, ttbdual,  tvsdual, 
    tposss, tpvdual, tvdual, kk, ttbintnumbs;
    #
    if giveupflag and Runtime()-st>giveupruntime then
      gvstopflag:=true;
      return();
    fi;
    leng:=Length(psol);
    #
    if leng=0 then beep(99131); fi;
    if leng=nn then 
      reachflag:=true;
      tg:=basisinv*psol;
      if not IsIntMat(tg) then beep(899112); fi;
      if TMTTmult(tg, GramL)<>GramL then beep(716552); fi;
      thetg:=List(tg);
      return();
    else 
      newsubcandidates:=[];
      ttbdual:=basisdual[leng+1];
      tvsdual:=candidatessdual[leng+1];
      ttbintnumbs:=List([1..leng], kk->newGram[leng+1][kk]);
      pos:=0;
      for tpvdual in tvsdual do 
        pos:=pos+1;
        ss:=3;
        for tvdual in [tpvdual, -tpvdual] do 
          ss:=ss-2; 
          # when tvdual=tpvdual then ss=1, when tvdual=-tpvdual then ss=-1, 
          if psol*tvdual=ttbintnumbs then 
            Add(newsubcandidates, [pos, ss]);
          fi;
        od;
      od;
      #
      if Length(newsubcandidates)<>fingerprints[leng+1] then 
        return();
      fi;
      #
      tvs:=candidatess[leng+1];
      for tposss in newsubcandidates do 
        tv:=tposss[2]*tvs[tposss[1]];
        Add(psol, tv);
        simpleextend(psol);
        if tv<>Remove(psol) then beep(812111); fi;
        if reachflag then return();; fi;
        if gvstopflag then return();; fi;
      od;
      if giveupflag then 
        if Runtime()-st>giveupruntime then 
          gvstopflag:=true;
        fi;
      fi;
    fi;
    return();
  end;
  #
  getorb:=function(aa, gens)
    local sorb,xx,tg,yy, orb;
    sorb:=[aa];
    for xx in sorb do 
      for tg in gens do 
        yy:=xx^tg;
        if not yy in sorb then Add(sorb, yy); fi;
      od;
    od;
    return(sorb);
  end;
  #
  #
  getorbsunion:=function(aas, gens)
    local orbunion, aa;
    orbunion:=[];
    for aa in aas do 
      if not aa in orbunion then 
        orbunion:=Union(orbunion, getorb(aa, gens));
      fi;
    od;
    return(orbunion);
  end;
  #
  subcandidates1:=subcandidatess[1];
  #
  inipos:=function(tv)
    local xx;
    for xx in tv do 
      if xx>0 then return(tv); elif xx<0 then return(-tv); fi;
    od;
    beep(727211);
  end;
  getpermcan1:=function(tg)
    local tvtgs, tv, tperm, ttv, tposs;
    tvtgs:=List(subcandidates1, tv->inipos(tv*tg));
    tperm:=PermList(List(tvtgs, ttv->Position(subcandidates1, ttv)));
    return(tperm);
  end;
  #
  st:=Runtime();
  #
  stabrecs:=[];
  #
  thelevel:=1;
  gengs:=[];
  genperms:=[];
  trueseeds:=[];
  falseseeds:=[];
  lengcan1:=Length(subcandidates1);
  doneposs:=[];
  pos:=0;
  for tv in subcandidates1 do 
    pos:=pos+1;
    if not pos in doneposs then 
      reachflag:=false;
      simpleextend([tv]);
      if gvstopflag then return(fail); fi;
      #
      if reachflag then 
        Add(gengs, thetg);
        tperm:=getpermcan1(thetg);
        Add(genperms, tperm);
        Add(trueseeds, pos);
        tlist:=getorbsunion(trueseeds, genperms);
        flist:=getorbsunion(falseseeds, genperms);
        doneposs:=Union(tlist, flist);
      else 
        Add(falseseeds, pos);
        neworb:=getorb(pos, genperms);
        doneposs:=Union(doneposs, neworb);
      fi;
      #
    fi;
  od;
  #
  trueposs:=getorbsunion(trueseeds, genperms);
  falseposs:=getorbsunion(falseseeds, genperms);
  if  Intersection(trueposs, falseposs)<>[]  then beep(18128181); fi;
  if Union(trueposs, falseposs)<>[1..lengcan1] then beep(18228181); fi;
  if trueposs=[] then beep(7766152); fi;
  #
  stabrec:=rec(
    level:=thelevel,
    candidates:=subcandidates1, 
    genperms:=genperms, 
    gengs:=gengs, 
    trueposs:=trueposs, 
    size:=Length(trueposs)
  );
  Add(stabrecs, stabrec);
  #
  # level 1 is done
  #
  for thelevel in [2..nn] do 
    #
    inipsol:=List([1..thelevel-1], kk->basis[kk]);
    subcandidates:=subcandidatess[thelevel];
    #
    getpermcan:=function(tg)
      local tvtgs, tv, tperm, ttv;
      tvtgs:=List(subcandidates, tv->tv*tg);
      tperm:=PermList(List(tvtgs, ttv->Position(subcandidates, ttv)));
      return(tperm);
    end;
    #
    gengs:=[];
    genperms:=[];
    trueseeds:=[];
    falseseeds:=[];
    lengcan:=Length(subcandidates);
    doneposs:=[];
    pos:=0;
    for tv in subcandidates do 
      pos:=pos+1;
      if not pos in doneposs then 
        reachflag:=false;
        simpleextend(CopyAdd(inipsol, tv));
        if gvstopflag then return(fail); fi;
        #
        if reachflag then 
          Add(gengs, thetg);
          tperm:=getpermcan(thetg);
          Add(genperms, tperm);
          Add(trueseeds, pos);
          tlist:=getorbsunion(trueseeds, genperms);
          flist:=getorbsunion(falseseeds, genperms);
          doneposs:=Union(tlist, flist);  
        else 
          Add(falseseeds, pos);
          neworb:=getorb(pos, genperms);
          doneposs:=Union(doneposs, neworb);
        fi;
        #
      fi;
    od;# for thelevel in [2..nn] do 
    #
    if gvstopflag then return(fail); fi;
    #
    trueposs:=getorbsunion(trueseeds, genperms);
    falseposs:=getorbsunion(falseseeds, genperms);
    if Intersection(trueposs, falseposs)<>[]  then beep(1822128181); fi;
    if Union(trueposs, falseposs)<>[1..lengcan] then beep(1822228181); fi;
    if trueposs=[] then beep(75516152); fi;
    stabrec:=rec(
      level:=thelevel,
      candidates:=subcandidates, 
      genperms:=genperms, 
      gengs:=gengs, 
      trueposs:=trueposs, 
      size:=Length(trueposs)
    );
    Add(stabrecs, stabrec);
  od;
  #
  totalgens:=[-IdentityMat(nn)];
  for stabrec in stabrecs do 
    Append(totalgens, stabrec.gengs);
  od;
  totalgens:=Set(totalgens);
  OGrec:=rec(
    Gram:=GramL,
    basis:=basis, 
    basisrec:=basisrec, 
    stabrecs:=stabrecs, 
    size:=2*Product(List(stabrecs, stabrec->stabrec.size)),
    totalgens:=totalgens
  );
  #
  return(OGrec);
  #
end;

CheckOGrecSize:=function(basisrec, OGrec)
  local vs, tg, tperm, gensperms, permsize;
  vs:=Union(basisrec.vss);
  vs:=Union(vs, -vs);
  gensperms:=[];
  for tg in OGrec.totalgens do
    tperm:=PermList(List(vs*tg, ttv->Position(vs, ttv)));
    Add(gensperms, tperm);
  od;
  permsize:=Size(Group(gensperms));
  if OGrec.size<>permsize then beep(9999111); fi;
  return(true);
end;


IsIsomBasisRecs:=function(basisrecA, basisrecB)
  #
  local beep, GramL1, GramL2, nn, 
  thevss1, thevss2, thevss1dual, basisrec1, basisrec2,
  fingerprints1, basis1,  basis1inv, totalflag, newGram1, newGram2, tt2vs, 
  thetg, extend, cleng, ii, kk, jj, wb1, wb2, revflag, thevss2dual, 
  webrints1, webrints2, basiswebrints1, candidatess2, bpos, tvs, 
  tcandidates, webrint2, tv, tvdual, candidatess2dual, vpos, ttvs;
  #
  beep:=function(beepnumb)
    localbeep("IsIsomBasisRec", beepnumb); Error();
  end;
  #
  if basisrecA.minflag and basisrecB.minflag then 
    if basisrecA.nrmsset<>basisrecB.nrmsset then return(false); fi;
    if basisrecA.nopss<>basisrecB.nopss then return(false); fi;
    revflag:=false;
    basisrec1:=basisrecA;
    basisrec2:=basisrecB;
    cleng:=Length(basisrecA.nrmsset);
  else 
    if Length(basisrecB.nrmsset)<=Length(basisrecA.nrmsset) and
       basisrecA.nrmsset{[1..Length(basisrecB.nrmsset)]}=basisrecB.nrmsset then 
      cleng:=Length(basisrecB.nrmsset);
      revflag:=true;
      basisrec1:=basisrecB;
      basisrec2:=basisrecA;
    elif Length(basisrecA.nrmsset)<=Length(basisrecB.nrmsset) and
         basisrecB.nrmsset{[1..Length(basisrecA.nrmsset)]}=basisrecA.nrmsset then 
      cleng:=Length(basisrecA.nrmsset);
      revflag:=false;
      basisrec1:=basisrecA;
      basisrec2:=basisrecB;
    else 
      return(false);
    fi;
    #
    for ii in [1..cleng] do 
      if basisrecA.nopss[ii]<>basisrecB.nopss[ii] then 
        return(false);
      fi;
    od;
  fi;
  #
  GramL1:=basisrec1.Gram;
  GramL2:=basisrec2.Gram;
  newGram1:=basisrec1.newGram;
  newGram2:=basisrec2.newGram;
  #
  nn:=Length(GramL1);
  if nn<>Length(GramL2) then return(false); fi;
  if DeterminantIntMat(GramL1)<>DeterminantIntMat(GramL2) then return(false); fi;
  #
  thevss1:=basisrec1.vss;
  thevss2:=List([1..cleng], kk->basisrec2.vss[kk]);
  thevss1dual:=List(thevss1, tvs->tvs*GramL1);
  thevss2dual:=List(thevss2, tvs->tvs*GramL2);
  #
  for jj in [1..cleng] do 
    for kk in [1..cleng] do
      wb1:=Collected(List(thevss1[kk], tv->Webrint(tv, GramL1, thevss1[jj])));
      wb2:=Collected(List(thevss2[kk], tv->Webrint(tv, GramL2, thevss2[jj])));
      if  wb1<>wb2 then return(false); fi;
    od;
  od;
  #
  basis1:=basisrec1.basis;
  basis1inv:=InverseMat(basis1);
  #
  basiswebrints1:=basisrec1.basiswebrintss;
  #
  candidatess2:=[];
  for bpos in [1..nn] do 
    vpos:=SinglePosition(basisrec1.nrmsset, basisrec1.basisnrms[bpos]);
    tcandidates:=[];
    for tv in thevss2[vpos] do 
      webrint2:=List(thevss2, tt2vs->Webrint(tv, GramL2, tt2vs));
      if webrint2=basiswebrints1[bpos] then 
        Add(tcandidates, tv);
      fi;
    od;
    #
    Add(candidatess2, tcandidates);
  od;
  candidatess2dual:=List(candidatess2, tvs->tvs*GramL2);
  #
  fingerprints1:=basisrec1.fingerprints;
  #
  totalflag:=false;
  thetg:=[];
  #
  extend:=function(psol)
    #
    local leng,tintnumbs, tvsdual, tvs, tvdual, pos,
    candidates, psoltvdual, tv, tsubcandidates;
    #
    leng:=Length(psol);
    if leng=nn then 
      totalflag:=true;
      thetg:=basis1inv*psol;
      if not IsIntMat(thetg) then beep(919191); fi;
      if TMTTmult(thetg, GramL2)<>GramL1 then beep(55112); fi;
      return();
    fi;
    #
    if leng=0 then
      tsubcandidates:=candidatess2[1];
    else 
      tsubcandidates:=[];
      tintnumbs:=List([1..leng], kk->newGram1[leng+1][kk]);
      tvsdual:=candidatess2dual[leng+1];
      tvs:=candidatess2[leng+1];
      pos:=0;
      for tvdual in tvsdual do 
        pos:=pos+1;
        psoltvdual:=psol*tvdual;
        if psoltvdual=tintnumbs then 
          Add(tsubcandidates, tvs[pos]);
        fi;
        #
        # this part cannot be changed to elif, 
        # because tintnumbs may be a zero vector
        #
        if psoltvdual=-tintnumbs then 
          Add(tsubcandidates, -tvs[pos]);
        fi;
      od;
    fi;
    #
    if Length(tsubcandidates)<>fingerprints1[leng+1] then return(); fi;
    #
    for tv in tsubcandidates do 
      Add(psol, tv);
      extend(psol);
      if Remove(psol)<>tv then beep(776611); fi;
      if totalflag then return(); fi;
    od;
    #
  end;
  #
  extend([]);
  #
  if totalflag then 
    if revflag then 
      thetg:=InverseMat(thetg);
    fi;
    if not IsIntMat(thetg) then beep(776611); fi;
    if TMTTmult(thetg, basisrecB.Gram)<>basisrecA.Gram then 
      beep(9137126); 
    fi;
    return(thetg);
  else 
    return(false);
  fi;
  #
end;


RepeatOGLatFromBasisRec:=function(GramL, basisrectrial, giveupsec)
  local counter, basisrec, OGrec;
  counter:=0;
  while true do
    counter:=counter+1;
    basisrec:=BasisRec(GramL,  basisrectrial+counter);
    OGrec:=OGLatFromBasisRec(basisrec, giveupsec+counter);
    if OGrec<>fail then 
      OGrec.basisrectrial:= basisrectrial+counter;
      OGrec.giveupsec:=giveupsec+counter;
      return(OGrec);
    else 
      Printn("___repeat", counter);
    fi;
  od;
end;

# OGLat(GramL):
# Returns a record with .size and .totalgens.
# The automorphism group is Group(result.totalgens).

# IsIsomLats(GramA, GramB):
# Returns false if the lattices are not isometric.
# Otherwise returns an integral unimodular matrix T satisfying
# TMTTmult(T, GramB) = GramA.



OGLat:=function(GramL)
  return(RepeatOGLatFromBasisRec(GramL, 3, 100));
end;

IsIsomLats:=function(GramLA, GramLB)
  local basisrecA, basisrecB;
  basisrecA:=BasisRec(GramLA, 10);
  basisrecB:=BasisRec(GramLB, 10);
  return(IsIsomBasisRecs(basisrecA, basisrecB));
end;



###########################################

#Read("SimpleOGLat.g");

#PleskenSouvignier

SimpleOGLat:=function(arg)
  local beep, GramL, montag, nn, LLLrec, basis, newGram, basisnrmslist,
  ii, maxnrm, svrec, nrmsset, pvposss, pvlists, aa, poss, pos, getpvlist,
  tbasis, tbasisdual, fingerprints, lllist, minll, inipos, 
  initb, OGrec, getfp,  tfps, tbs, tb, minfp, newtb,
  basisinv, basisdual, basisnrms, reachflag, thetg, simpleextend,
  getorb, candidates1, getpermcan1, stabrecs, stabrec, doneposs,
  gengs, genperms, trueposs, neworb, tv, candidates, getpermcan, tperm,
  inipsol, kk, GramLinv, pvduallists, getpvduallist, tpvlist,
  getorbsunion, falseseeds, trueseeds, falseposs, thelevel, inipsoldual,
  tintnumbs, tpv, jj, intnumbss, totalgens, lengcan, lengcan1, tblist, fblist,
  getorbblist, getorbsunionblist, colnrms;
  #
  beep:=function(beepnumb)
    localbeep("SimpleOGLat", beepnumb); Error();
  end;
  #
  GramL:=arg[1]; 
  if Length(arg)>1 then montag:=arg[2]; else montag:=0; fi;
  nn:=Length(GramL);
  GramLinv:=InverseMat(GramL);
  LLLrec:=LLLReducedGramMat(GramL);
  basis:=LLLrec.transformation;
  newGram:=LLLrec.remainder;
  if newGram<>TMTTmult(basis, GramL) then beep(11211); fi;
  basisnrmslist:=List([1..nn], ii->newGram[ii][ii]);
  maxnrm:=Maximum(basisnrmslist);
  svrec:=ShortestVectors(GramL, maxnrm);
  nrmsset:=Set(svrec.norms);
  colnrms:=Collected(svrec.norms);
  pvposss:=List(nrmsset, aa->Positions(svrec.norms, aa));
  #
  inipos:=function(tv)
    local xx;
    for xx in tv do 
      if xx>0 then return(tv); 
      elif xx<0 then return(-tv);
      fi;
    od;
    beep(727211);
  end;
  #
  #pvlists:=List(pvposss, poss->List(poss, pos->inipos(svrec.vectors[pos])));
  pvlists:=List(pvposss, poss->Set(poss, pos->inipos(svrec.vectors[pos])));
  pvduallists:=List(pvlists, tpvlist->tpvlist*GramL);
  #
  getpvlist:=function(nrm)
    return(pvlists[SinglePosition(nrmsset, nrm)]);
  end;
  #
  getpvduallist:=function(nrm)
    return(pvduallists[SinglePosition(nrmsset, nrm)]);
  end;
  #
  tbasis:=[];
  tbasisdual:=[];
  fingerprints:=[];
  #
  lllist:=List(basisnrmslist, aa->Length(getpvlist(aa)));
  minll:=Minimum(lllist);
  initb:=basis[Position(lllist,  minll)];
  Add(fingerprints, minll);
  Add(tbasis, initb);
  Add(tbasisdual, initb*GramL);
  #
  getfp:=function(tb)
    local tintnumbs, ttintnumbs, tpvlist, fp, tpv, tv;
    tintnumbs:=tbasisdual*tb;
    tpvlist:=getpvlist(tb*GramL*tb);
    fp:=0;
    for tpv in tpvlist do 
      ttintnumbs:=tbasisdual*tpv;
      if ttintnumbs=tintnumbs  then fp:=fp+1; fi;
      if ttintnumbs=-tintnumbs  then fp:=fp+1; fi;
    od;
    return(fp);
  end;
  #
  for jj in [2..nn] do 
    tfps:=[];
    tbs:=[];
    for tb in basis do 
      if tb in tbasis then continue; fi;
      Add(tfps, getfp(tb));
      Add(tbs, tb);
    od;
    minfp:=Minimum(tfps);
    Add(fingerprints,  minfp);
    newtb:=tbs[Position(tfps, minfp)];
    Add(tbasis, newtb);
    Add(tbasisdual, newtb*GramL);
    if montag>=2 then 
      Printn("___making fingerprints", jj);
    fi;
  od;
  if montag>=1 then 
    Printn("___fingerprints", fingerprints);
  fi;
  #
  if not IsEqualSet(basis, tbasis) then beep(818121); fi;
  basis:=ShallowCopy(tbasis);
  basisinv:=InverseMat(basis);
  basisdual:=basis*GramL;
  basisnrms:=List(basis, tb->tb*GramL*tb);
  newGram:=TMTTmult(basis, GramL);
  intnumbss:=List([1..nn-1], jj->List([1..jj], kk->newGram[jj+1][kk]));
  #
  reachflag:=false;
  thetg:=[];
  #
  simpleextend:=function(psol)
    local leng, tg, candidates,  tpvlist, intnumbs,
    tpv, tv, tintnumbs, pos, ss, ttbdual,  tpvduallist, 
    tposss, tpvdual, tvdual;
    #psoldual; #bdg
    #
    #
    leng:=Length(psol);
    #
    if leng=0 then beep(99131); fi;
    if leng=nn then 
      reachflag:=true;
      tg:=basisinv*psol;
      if not IsIntMat(tg) then beep(899112); fi;
      if TMTTmult(tg, GramL)<>GramL then beep(716552); fi;
      thetg:=List(tg);
      return();
    else 
      candidates:=[];
      ttbdual:=basisdual[leng+1];
      tpvduallist:=getpvduallist(basisnrms[leng+1]);
      tintnumbs:=intnumbss[leng];
      pos:=0;
      for tpvdual in tpvduallist do 
        pos:=pos+1;
        ss:=3;
        for tvdual in [tpvdual,-tpvdual] do 
          ss:=ss-2;
          if tintnumbs=psol*tvdual then 
            Add(candidates, [pos, ss]);
          fi;
        od;
      od;
      if Length(candidates)<>fingerprints[leng+1] then 
        return();
      fi;
      tpvlist:=getpvlist(basisnrms[leng+1]);
      for tposss in candidates do 
        tv:=tposss[2]*tpvlist[tposss[1]];
        Add(psol, tv);
        simpleextend(psol);
        if tv<>Remove(psol) then beep(812111); fi;
        if reachflag then return();; fi;
      od;
    fi;
    return();
  end;
  #
  getorbblist:=function(aa, gens, theleng)
    local sorb,xx,tg,yy, orb;
    sorb:=[aa];
    for xx in sorb do 
      for tg in gens do 
        yy:=xx^tg;
        if not yy in sorb then Add(sorb, yy); fi;
        # Naive breadth-first orbit computation
      od;
    od;
    orb:=BlistList([1..theleng], sorb);
    return(orb);
  end;
  #
  # getorbblist:=function(aa, gens, theleng)
  #   local orb,task;
  #   orb:=BlistList([1..theleng], []);
  #   task:=function(xx)
  #     local tg;
  #     if not orb[xx] then 
  #       orb[xx]:=true;
  #       for tg in gens do 
  #         task(xx^tg);
  #       od;
  #     fi;
  #   end;
  #   task(aa);
  #   return(orb);
  # end;
  #
  # getorbsunion:=function(aas, gens)
  #   local orbunion, aa;
  #   orbunion:=[];
  #   for aa in aas do 
  #     if not orbunion[aa] then 
  #       orbunion:=Union(orbunion, getorb(aa, gens));
  #     fi;
  #   od;
  #   return(orbunion);
  # end;
  #
  getorbsunionblist:=function(aas, gens, theleng)
    local orbunion, aa;
    orbunion:=BlistList([1..theleng], []);;
    for aa in aas do 
      if not aa in orbunion then 
        UniteBlist(orbunion, getorbblist(aa, gens, theleng));
      fi;
    od;
    return(orbunion);
  end;
  #
  candidates1:=getpvlist(basisnrms[1]);
  #
  getpermcan1:=function(tg)
    local tvtgs, tv, tperm, ttv, tposs;
    tvtgs:=List(candidates1, tv->inipos(tv*tg));
    #tperm:=PermList(List(tvtgs, ttv->SinglePosition(candidates1, ttv)));
    #tposs:=List(tvtgs, ttv->Position(candidates1, ttv));
    #if fail in tposs then beep(442511); fi;
    tperm:=PermList(List(tvtgs, ttv->Position(candidates1, ttv)));
    return(tperm);
  end;
  #
  stabrecs:=[];
  #
  thelevel:=1;
  gengs:=[];
  genperms:=[];
  trueseeds:=[];
  falseseeds:=[];
  lengcan1:=Length(candidates1);
  doneposs:=BlistList([1..lengcan1], []);
  pos:=0;
  for tv in candidates1 do 
    pos:=pos+1;
    if not doneposs[pos] then 
      reachflag:=false;
      simpleextend([tv]);
      #
      if reachflag then 
        Add(gengs, thetg);
        tperm:=getpermcan1(thetg);
        Add(genperms, tperm);
        Add(trueseeds, pos);
        tblist:=getorbsunionblist(trueseeds, genperms, lengcan1);
        fblist:=getorbsunionblist(falseseeds, genperms, lengcan1);
        doneposs:=tblist;
        UniteBlist(doneposs, fblist);
      else 
        Add(falseseeds, pos);
        neworb:=getorbblist(pos, genperms, lengcan1);
        UniteBlist(doneposs, neworb);
      fi;
      #
    fi;
  od;
  #
  trueposs:=getorbsunionblist(trueseeds, genperms, lengcan1);
  falseposs:=getorbsunionblist(falseseeds, genperms, lengcan1);
  if true in IntersectionBlist(trueposs, falseposs)then beep(18128181); fi;
  if false in UnionBlist(trueposs, falseposs) then beep(18228181); fi;
  trueposs:=Set(Positions(trueposs, true));
  #
  stabrec:=rec(
    level:=thelevel,
    candidates:=candidates1, 
    genperms:=genperms, 
    gengs:=gengs, 
    trueposs:=trueposs, 
    order:=Length(trueposs)
  );
  Add(stabrecs, stabrec);
  if montag>=1 then 
    Printn("level", thelevel, stabrec.order);
  fi;
  #
  #
  for thelevel in [2..nn] do 
    #
    #
    inipsol:=List([1..thelevel-1], kk->basis[kk]);
    inipsoldual:=inipsol*GramL;
    tintnumbs:=intnumbss[thelevel-1];
    tpvlist:=getpvlist(basisnrms[thelevel]);
    candidates:=[];
    for tpv in tpvlist do 
      for tv in [tpv, -tpv] do 
        if tintnumbs=inipsoldual*tv then 
          Add(candidates, tv);
        fi;
      od;
    od;
    #
    getpermcan:=function(tg)
      local tvtgs, tv, tperm, ttv;
      tvtgs:=List(candidates, tv->tv*tg);
      #tperm:=PermList(List(tvtgs, ttv->SinglePosition(candidates, ttv)));
      tperm:=PermList(List(tvtgs, ttv->Position(candidates, ttv)));
      return(tperm);
    end;
    #
    gengs:=[];
    genperms:=[];
    trueseeds:=[];
    falseseeds:=[];
    lengcan:=Length(candidates);
    doneposs:=BlistList([1..lengcan], []);
    pos:=0;
    for tv in candidates do 
      pos:=pos+1;
      if not doneposs[pos] then 
        reachflag:=false;
        simpleextend(CopyAdd(inipsol, tv));
        #
        if reachflag then 
          Add(gengs, thetg);
          tperm:=getpermcan(thetg);
          Add(genperms, tperm);
          Add(trueseeds, pos);
          tblist:=getorbsunionblist(trueseeds, genperms, lengcan);
          fblist:=getorbsunionblist(falseseeds, genperms, lengcan);
          doneposs:=tblist;
          UniteBlist(doneposs, fblist);
        else 
          Add(falseseeds, pos);
          neworb:=getorbblist(pos, genperms, lengcan);
          UniteBlist(doneposs, neworb);
        fi;
        #
      fi;
    od;
    #
    trueposs:=getorbsunionblist(trueseeds, genperms, lengcan);
    falseposs:=getorbsunionblist(falseseeds, genperms, lengcan);
    if true in IntersectionBlist(trueposs, falseposs)then beep(8128181); fi;
    if false in UnionBlist(trueposs, falseposs) then beep(8228181); fi;
    trueposs:=Set(Positions(trueposs, true));
    stabrec:=rec(
      level:=thelevel,
      candidates:=candidates, 
      genperms:=genperms, 
      gengs:=gengs, 
      trueposs:=trueposs, 
      order:=Length(trueposs)
    );
    Add(stabrecs, stabrec);
    if montag>=1 then 
      Printn("level", thelevel, stabrec.order);
    fi;
  od;
  #
  totalgens:=[-IdentityMat(nn)];
  for stabrec in stabrecs do 
    Append(totalgens, stabrec.gengs);
  od;
  totalgens:=Set(totalgens);
  OGrec:=rec(
    Gram:=GramL,
    det:=DeterminantIntMat(GramL),
    colnorms:=colnrms, 
    basis:=basis, 
    stabrecs:=stabrecs, 
    order:=2*Product(List(stabrecs, stabrec->stabrec.order)),
    totalgens:=totalgens
  );
  #
  return(OGrec);
  #
end;




GensSize:=function(GramL, gens)
  local pvs, vs, permgens, tnrm, tb, tpg, tg, tbasis, tv, size;
  #
  if Set(gens, tg->TMTTmult(tg, GramL)=GramL)<>[true] then 
    beep(666555); 
  fi;
  tbasis:=LLLReducedGramMat(GramL).transformation;
  tnrm:=Maximum(List(tbasis, tb->tb*GramL*tb));
  pvs:=ShortestVectors(GramL, tnrm).vectors;
  vs:=Union(pvs, -pvs);
  permgens:=[];
  for tg in gens do 
    tpg:=PermList(List(vs*tg, tv->Position(vs, tv)));
    Add(permgens, tpg);
  od;
  size:=Size(Group(permgens));
  #
  return(size);
  #
end;


LLLNormByRandU:=function(Gram)
  local nn, nrms, rU, tbasis, tnrm, tGram;
  nn:=Length(Gram);
  nrms:=[];
  while Length(nrms) <10 do
    rU:=RandomUnimodMat(nn);
    tGram:=TMTTmult(rU, Gram);
    tbasis:=LLLReducedGramMat(tGram).transformation;
    tnrm:=Maximum(List(tbasis, tb->tb*tGram*tb));
    Add(nrms, tnrm);
  od;
  return(Minimum(nrms));
end;

NewGensSize:=function(GramL, gens)
  local pvs, vs, permgens, tnrm, tb, tpg, tg, tbasis, tv, size;
  #
  if Set(gens, tg->TMTTmult(tg, GramL)=GramL)<>[true] then 
    beep(666555); 
  fi;
  tnrm:=LLLNormByRandU(GramL);
  pvs:=ShortestVectors(GramL, tnrm).vectors;
  vs:=Union(pvs, -pvs);
  permgens:=[];
  for tg in gens do 
    tpg:=PermList(List(vs*tg, tv->Position(vs, tv)));
    Add(permgens, tpg);
  od;
  size:=Size(Group(permgens));
  #
  return(size);
  #
end;


#
# Made by ChatGPT on 2026/09/28 
#
# Requires GAP's LLLReducedGramMat and ShortestVectors.
# OGrecL1 is the existing, unmodified result of SimpleOGLat(GramL1).
# Positive definite integral Gram matrices of positive rank are assumed.
# Returns fail iff the lattices are not isomorphic; otherwise returns T with
#   T*GramL2*TransposedMat(T) = OGrecL1.Gram and Det(T) = +/-1.
# In row coordinates the isomorphism L1 -> L2 is v -> v*T.
# Neither SimpleOGLat nor an automorphism computation for L2 is called.

SimpleIsomLats:=function(OGrecL1, GramL2)
  #
  local GramL1, n, basis1, basis1inv, sourceGram, norms,
  red2, basis2, targetGram, sv, shells, duals, fingerprints,
  positive, search, result, first, getshell, beep, ss, vv;
  #
  beep:=function(beepnumb)
    localbeep("SimpleIsomLats", beepnumb); Error();
  end;
  #
  GramL1:=OGrecL1.Gram;
  n:=Length(GramL1);
  if n<>Length(GramL2) then return fail; fi;
  if n=0 then Error("Positive rank required"); fi;
  if OGrecL1.det<>DeterminantIntMat(GramL2) then
    return fail;
  fi;
  if GramL1=GramL2 then return IdentityMat(n); fi;
  #
  basis1:=OGrecL1.basis;
  basis1inv:=InverseMat(basis1);
  sourceGram:=basis1*GramL1*TransposedMat(basis1);
  norms:=List([1..n],i->sourceGram[i][i]);
  # Level 1 candidates are modulo signs; later levels are signed.
  # These lengths recover the fingerprints from the existing SimpleOGLat record.
  fingerprints:=List(OGrecL1.stabrecs,r->Length(r.candidates));
  # Enumerate in reduced target coordinates; convert back on return.
  #
  red2:=LLLReducedGramMat(GramL2);
  basis2:=red2.transformation;
  targetGram:=red2.remainder;
  sv:=ShortestVectors(targetGram, Maximum(norms));
  if OGrecL1.colnorms<>Collected(sv.norms) then return fail; fi;
  #
  positive:=function(v)
    local x;
    for x in v do
      if x>0 then return v; elif x<0 then return -v; fi;
    od;
    Error("Unexpected zero short vector");
  end;
  #
  getshell:=function(a)
    local poss, pos;
    # positions of norm a vectors 
    poss:=Filtered([1..Length(sv.norms)], pos->sv.norms[pos]=a);
    # sign-normalize and make duplicate-free 
    return Set(poss, pos->positive(sv.vectors[pos]));
  end;
  #
  shells:=List(norms, getshell);
  #
  if Length(shells[1])<>fingerprints[1] then return fail; fi;
  if ForAny(shells, ss->Length(ss)=0) then return fail; fi;
  duals:=List(shells, ss->List(ss, vv->vv*targetGram));
  result:=fail;
  #
  search:=function(rows)
    local k, wanted, candidates, p, products, v, T;
    k:=Length(rows)+1;
    if k>n then
      T:=basis1inv*rows*basis2;
      if not IsIntMat(T) then beep(76154361); fi;
      if AbsInt(DeterminantIntMat(T))<>1 then beep(7615761); fi;
      if TMTTmult(T, GramL2)<>GramL1 then beep(763331); fi;
      result:=T;
      return;
    fi;
    wanted:=sourceGram[k]{[1..k-1]};
    candidates:=[];
    for p in [1..Length(shells[k])] do
      products:=rows*duals[k][p];
      if products=wanted then Add(candidates, shells[k][p]); fi;
      if -products=wanted then Add(candidates, -shells[k][p]); fi;
    od;
    if Length(candidates)<>fingerprints[k] then return; fi;
    for v in candidates do
      Add(rows,v);
      search(rows);
      if v<>Remove(rows) then beep(4714761); fi;
      if result<>fail then return; fi;
    od;
  end;
  #
  # Any isometry can be multiplied by -1; one sign per first image suffices.
  for first in shells[1] do
    search([first]);
    if result<>fail then return result; fi;
  od;
  return fail;
end;



#####