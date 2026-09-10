#Read("NewOGLat.g");


VssRec:=function(GramL, nrmcandidate)
  #
  # If nrmcandidate is the minimal norm such that
  # the vevtors v with v^2 <=nrmcandidate generate L,
  # then this returms vssrec with vssrec.minflag:=true;
  # elif nrmcandidate is too large (not minimal),  then 
  # this returms vssrec with vssrec.minflag:=false;
  # elif nrmcandidate is too small,  then it beeps.
  #
  local nn, svsrec, nrmsset, vects, vss, tvs, tnrm, 
  tpos, doesgenerate, vssrec, poss, vs, minflag,
  minreqnrm; #minimal required norm 
  #
  nn:=Length(GramL);
  svsrec:=ShortestVectors(GramL, nrmcandidate);
  nrmsset:=Set(svsrec.norms);
  if Maximum(nrmsset)<nrmcandidate then return(false); fi;
  vects:=svsrec.vectors;
  vss:=[];
  tvs:=[];
  minflag:=true; 
  minreqnrm:=nrmcandidate;
  for tnrm  in nrmsset do 
    poss:=Positions(svsrec.norms, tnrm);
    vs:=List(poss, tpos->vects[tpos]);
    Add(vss, vs);
    Append(tvs, vs);
    doesgenerate:=(Rank(tvs)=nn and CokerTorsion(tvs)=[]);
    if doesgenerate and tnrm<nrmcandidate then 
      minflag:=false; 
      minreqnrm:=tnrm;
    fi;
    if tnrm=nrmcandidate and (not doesgenerate) then 
      Printn("does not generate!");
      beep(442211); 
    fi;
  od;
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


BasisRec:=function(GramL, trialnumb)
  #
  local nn, norm, tU, tt,GramL2,  LLLrec, basis, ii,  
  basisnrms,fflag,  maxnrm,  tbasisnrms, 
  newGram, getvs, trec, jj, kk, tintnumbs,
  tbasisdual, tv, tvs,  tbasis, thebasis,
  nrmcandidate, ttintnumbs, maxbasisrrms,  tmaxbasisnrms,
  counter, minreqnrm, thevss,  basiswebrintss,
  candidatess, bpos, tcandidates, inipos, subcandidatess, tsubcandidates,
  tbintnumbs, fingerprints;
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
  newGram:=TMTTmult(thebasis, GramL);
  #
  trec.basisnrms:=basisnrms;
  trec.basis:=thebasis;
  trec.newGram:=newGram;
  #
  thevss:=trec.vss;
  #
  #
  basiswebrintss:=List(thebasis, tb->Webrints(tb, GramL, thevss));
  #
  candidatess:=[];
  for bpos in [1..nn] do 
    tvs:=thevss[SinglePosition(trec.nrmsset,trec.basisnrms[bpos])];
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
  inipos:=function(tv)
    local xx;
    for xx in tv do 
      if xx>0 then return(tv); elif xx<0 then return(-tv); fi;
    od;
    beep(727211);
  end;
  candidatess[1]:=List(candidatess[1], inipos);
  # to make the computation in getpermcan1 easy
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
  trec.basiswebrintss:=basiswebrintss;
  trec.candidatess:=candidatess;
  trec.subcandidatess:=subcandidatess;
  trec.fingerprints:=fingerprints;
  #
  return(trec);
end;


NewOGLat:=function(arg) # arg is (basisrec) or (basisrec, giveupsec)
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
    localbeep("NewOGLat", beepnumb); Error();
  end;
  #
  if Length(arg)=1 then 
    basisrec:=arg[1];
    giveupflag:=false;
    giveupruntime:=infinity;
  elif Length(arg)=2 then 
    basisrec:=arg[1];
    giveupfla


g:=true;
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
    if IsSubset(basisrecA.nrmsset, basisrecB.nrmsset) then 
      cleng:=Length(basisrecB.nrmsset);
      revflag:=true;
      basisrec1:=basisrecB;
      basisrec2:=basisrecA;
    elif IsSubset(basisrecB.nrmsset, basisrecA.nrmsset) then 
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
    tcandidates:=candidatess2[leng+1];
    if leng=0 then 
      tsubcandidates:= tcandidates;
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


RepeatNewOGLat:=function(GramL, basisrectrial, giveupsec)
  local counter, basisrec, OGrec;
  counter:=0;
  while true do
    counter:=counter+1;
    basisrec:=BasisRec(GramL,  basisrectrial+counter);
    OGrec:=NewOGLat(basisrec, giveupsec+counter);
    if OGrec<>fail then 
      OGrec.basisrectrial:= basisrectrial+counter;
      OGrec.giveupsec:=giveupsec+counter;
      return(OGrec);
    else 
      Printn("___repeat", counter);
    fi;
  od;
end;

#####