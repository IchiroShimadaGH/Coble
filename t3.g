#Read("NewOGLat.g");


VssRec:=function(GramL, nrmcandidate)
  #
  # If nrmcandidate is the minimal norm such that
  # the vevtors v with v^2 <=nrmcandidate generate L,
  # then this returms vssrec;
  # elif nrmcandidate is too large (not minimal)  then it return (false);
  # elif nrmcandidate is too small,  then it beeps.
  #
  local nn, svsrec, nrmsset, vects, vss, tvs, tnrm, 
  tpos, doesgenerate, vssrec, poss, vs, dbg1;
  #
  nn:=Length(GramL);
  svsrec:=ShortestVectors(GramL, nrmcandidate);
  nrmsset:=Set(svsrec.norms);
  if Maximum(nrmsset)<nrmcandidate then return(false); fi;
  vects:=svsrec.vectors;
  vss:=[];
  tvs:=[];
  for tnrm in nrmsset do 
    poss:=Positions(svsrec.norms, tnrm);
    vs:=List(poss, tpos->vects[tpos]);
    Add(vss, vs);
    Append(tvs, vs);
    doesgenerate:=(Rank(tvs)=nn and CokerTorsion(tvs)=[]);
    if doesgenerate and tnrm<nrmcandidate then return(false); fi;
    if tnrm=nrmcandidate and (not doesgenerate) then 
      Printn("does not generate!");
      beep(442211); 
    fi;
  od;
  #
  dbg1:=List(vss, tvs->Set(tvs, tv->tv*GramL*tv));
  if dbg1<>List(nrmsset, tnrm->[tnrm] ) then 
    beep(2222266611111); #bdg
  fi;
  #
  vssrec:=rec(
    Gram:=List(GramL), 
    maxnrm:=nrmcandidate,
    nrmsset:=nrmsset, 
    nopss:=List(vss, Length),
    vss:=vss
  );
  return(vssrec);
end;





BasisRec:=function(GramL, trialnumb)
  #
  local nn, norm, tU, tt,GramL2,  LLLrec, basis, ii,  
  basisnrms,fflag,  maxnrm,  tbasisnrms, 
  fingerprints, newGram, getvs, trec, jj, kk, tintnumbs,
  tbasisdual, tv, tvs, fp, tbasis, thebasis,
  nrmcandidate, ttintnumbs, maxbasisrrms,  tmaxbasisnrms;
  #
  nn:=Length(GramL);
  if SignatureQ(GramL)<>[nn, nn, 0] then beep(665522); fi;
  fflag:=false;
  basisnrms:=[];
  maxbasisrrms:=infinity;
  while not fflag do
    for tt in [1..trialnumb] do
      tU:=RandomUnimodMat(nn);
      GramL2:=TMTTmult(tU, GramL);
      LLLrec:=LLLReducedGramMat(GramL2);
      tbasis:=LLLrec.transformation;
      #newGram:=LLLrec.remainder;
      #if newGram<>TMTTmult(basis, GramL2) then beep(11211); fi;
      tbasisnrms:=List([1..nn], ii->LLLrec.remainder[ii][ii]);
      tmaxbasisnrms:=Maximum(tbasisnrms);
      if basisnrms=[] or 
        maxbasisrrms> tmaxbasisnrms or 
       (maxbasisrrms= tmaxbasisnrms and basisnrms>tbasisnrms) then
        thebasis:=tbasis*tU;
        basisnrms:=List(tbasisnrms);
        maxbasisrrms:= tmaxbasisnrms;
      fi;
    od;
    nrmcandidate:=maxbasisrrms;
    trec:=VssRec(GramL, nrmcandidate);
    if trec<>false then fflag:=true; fi;
  od;
  #
  getvs:=function(tnrm)
    return(trec.vss[SinglePosition(trec.nrmsset, tnrm)]);
  end;
  #
  fingerprints:=[Length(getvs(basisnrms[1]))];
  #
  newGram:=TMTTmult(thebasis, GramL);
  if List([1..nn], kk->newGram[kk][kk])<>basisnrms then beep(52213); fi;
  #
  for jj in [2..nn] do 
    tintnumbs:=List([1..jj-1], kk->newGram[jj][kk]);
    tvs:=getvs(basisnrms[jj]);
    tbasisdual:=List([1..jj-1], kk->thebasis[kk])*GramL;
    fp:=0;
    for tv in tvs do 
      ttintnumbs:=tbasisdual*tv;
      if ttintnumbs=tintnumbs  then fp:=fp+1; fi;
      #
      # this part cannot be changed to elif, because tintnumbs may be a zero vector
      #
      if ttintnumbs=-tintnumbs  then fp:=fp+1; fi;
    od;
    Add(fingerprints, fp);
  od;
  #
  #
  trec.basisnrms:=basisnrms;
  trec.basis:=thebasis;
  trec.newGram:=newGram;
  trec.fingerprints:=fingerprints;
  #
  return(trec);
end;


NewOGLat:=function(basisrec)
  local beep,  GramL, nn,  basis, newGram, basisnrms,
  ii, maxnrm, svrec, nrmsset, pvposss, vss, aa, poss, pos, getvs,
  tbasis, tbasisdual, fingerprints,inipos, thevss, thevssdual, 
  initb, OGrec, getfp,  tfps, tbs, tb, minfp, newtb,
  basisinv, basisdual,  reachflag, thetg, simpleextend,
  getorb, candidates1, getpermcan1, stabrecs, stabrec, doneposs,
  gengs, genperms, trueposs, neworb, tv, candidates, getpermcan, tperm,
  inipsol, kk,  tvduals, getvsdual, tvs,
  getorbsunion, falseseeds, trueseeds, falseposs, thelevel, inipsoldual,
  tintnumbs, tpv, jj,  totalgens, lengcan, lengcan1, tblist, fblist,
  getorbblist, getorbsunionblist, fpflag, dbg1, pos1;
  #
  beep:=function(beepnumb)
    localbeep("NewOGLat", beepnumb); Error();
  end;
  #
  GramL:=basisrec.Gram;
  nn:=Length(GramL);
  basis:=basisrec.basis;
  fingerprints:=basisrec.fingerprints;
  newGram:=basisrec.newGram;
  if TMTTmult(basis, GramL)<>newGram then beep(391919); fi;
  nrmsset:=basisrec.nrmsset;
  basisnrms:=basisrec.basisnrms;
  basisinv:=InverseMat(basis);
  basisdual:=basis*GramL;
  #
  inipos:=function(tv)
    local xx;
    for xx in tv do 
      if xx>0 then return(tv); elif xx<0 then return(-tv); fi;
    od;
    beep(727211);
  end;
  #
  thevss:=basisrec.vss;
  pos1:=SinglePosition(nrmsset, basisnrms[1]);
  thevss[pos1]:=List(thevss[pos1], inipos);
  # to make the computation in getpermcan1 easy 
  #
  thevssdual:=List(thevss, tvs->tvs*GramL);
  dbg1:=List(thevss, tvs->Set(tvs, tv->tv*GramL*tv));
  if dbg1<>List(nrmsset, tnrm->[tnrm] ) then 
    beep(66611111); #bdg
  fi;
  #
  getvs:=function(nrm)
    return(thevss[SinglePosition(nrmsset, nrm)]);
  end;
  #
  getvsdual:=function(nrm)
    return(thevssdual[SinglePosition(nrmsset, nrm)]);
  end;
  #
  #
  reachflag:=false;
  thetg:=[];
  #
  simpleextend:=function(psol)
    local leng, tg, candidates,  tvs, 
    tpv, tv, tintnumbs, pos, ss, ttbdual,  tvsdual, 
    tposss, tpvdual, tvdual, kk;
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
      tvsdual:=getvsdual(basisnrms[leng+1]);
      tintnumbs:=List([1..leng], kk->newGram[leng+1][kk]);
      pos:=0;
      for tpvdual in tvsdual do 
        pos:=pos+1;
        ss:=3;
        for tvdual in [tpvdual, -tpvdual] do 
          ss:=ss-2; 
          # when tvdual=tpvdual then ss=1, when tvdual=-tpvdual then ss=-1, 
          if tintnumbs=psol*tvdual then 
            Add(candidates, [pos, ss]);
          fi;
        od;
      od;
      if thelevel=1 then #bdg
        Printn("in simpleextend", leng+1, Length(candidates));
      fi;
      tvs:=getvs(basisnrms[leng+1]);
      for tposss in candidates do 
        tv:=tposss[2]*tvs[tposss[1]];
        Add(psol, tv);
        simpleextend(psol);
        if tv<>Remove(psol) then beep(812111); fi;
        if reachflag then return();; fi;
      od;
    fi;
    return();
  end;
  #
  getorbblist:=function(aa, gens)
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
  getorbsunionblist:=function(aas, gens)
    local orbunion, aa;
    orbunion:=[];
    for aa in aas do 
      if not aa in orbunion then 
        orbunion:=Union(orbunion, getorbblist(aa, gens));
      fi;
    od;
    return(orbunion);
  end;
  #
  candidates1:=getvs(basisnrms[1]);
  #
  getpermcan1:=function(tg)
    local tvtgs, tv, tperm, ttv, tposs;
    tvtgs:=List(candidates1, tv->inipos(tv*tg));
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
  doneposs:=[];
  pos:=0;
  for tv in candidates1 do 
    pos:=pos+1;
    if not pos in doneposs then 
      reachflag:=false;
      simpleextend([tv]);
      Printn(pos, reachflag);
      #
      if reachflag then 
        Add(gengs, thetg);
        tperm:=getpermcan1(thetg);
        Add(genperms, tperm);
        Add(trueseeds, pos);
        tblist:=getorbsunionblist(trueseeds, genperms);
        fblist:=getorbsunionblist(falseseeds, genperms);
        doneposs:=Union(tblist, fblist);
      else 
        Add(falseseeds, pos);
        neworb:=getorbblist(pos, genperms);
        doneposs:=Union(doneposs, neworb);
      fi;
      #
    fi;
  od;
  #
  trueposs:=getorbsunionblist(trueseeds, genperms);
  falseposs:=getorbsunionblist(falseseeds, genperms);
  if  Intersection(trueposs, falseposs)<>[]  then beep(18128181); fi;
  if Union(trueposs, falseposs)<>[1..lengcan1] then beep(18228181); fi;
  if trueposs=[] then beep(7766152); fi;
  #
  stabrec:=rec(
    level:=thelevel,
    candidates:=candidates1, 
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
    inipsoldual:=inipsol*GramL;
    tintnumbs:=List([1..thelevel-1], kk->newGram[thelevel][kk]);
    tvs:=getvs(basisnrms[thelevel]);
    candidates:=[];
    for tpv in tvs do 
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
      tperm:=PermList(List(tvtgs, ttv->Position(candidates, ttv)));
      return(tperm);
    end;
    #
    gengs:=[];
    genperms:=[];
    trueseeds:=[];
    falseseeds:=[];
    lengcan:=Length(candidates);
    doneposs:=[];
    pos:=0;
    for tv in candidates do 
      pos:=pos+1;
      if not pos in doneposs then 
        reachflag:=false;
        simpleextend(CopyAdd(inipsol, tv));
        #
        if reachflag then 
          Add(gengs, thetg);
          tperm:=getpermcan(thetg);
          Add(genperms, tperm);
          Add(trueseeds, pos);
          tblist:=getorbsunionblist(trueseeds, genperms);
          fblist:=getorbsunionblist(falseseeds, genperms);
          doneposs:=Union(tblist, fblist);  
        else 
          Add(falseseeds, pos);
          neworb:=getorbblist(pos, genperms);
          doneposs:=Union(doneposs, neworb);
        fi;
        #
      fi;
    od;
    #
    trueposs:=getorbsunionblist(trueseeds, genperms);
    falseposs:=getorbsunionblist(falseseeds, genperms);
    if Intersection(trueposs, falseposs)<>[]  then beep(1822128181); fi;
    if Union(trueposs, falseposs)<>[1..lengcan] then beep(1822228181); fi;
    if trueposs=[] then beep(75516152); fi;
    stabrec:=rec(
      level:=thelevel,
      candidates:=candidates, 
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


FindIsomBasisRecs:=function(basisrec1, basisrec2)
  #
  local GramL1, GramL2, nn, vss2, vss2dual, getvs, getvsdual, basis1, fingerprints1, basis1inv, totalflag, newGram1, 
  thetg, extend;
  #
  if basisrec1.nrmsset<>basisrec2.nrmsset then return(false); fi;
  if basisrec1.nopss<>basisrec2.nopss then return(false); fi;
  #
  GramL1:=basisrec1.Gram;
  GramL2:=basisrec2.Gram;
  nn:=Length(GramL1);
  if nn<>Length(GramL2) then return(false); fi;
  if DeterminantIntMat(GramL1)<>DeterminantIntMat(GramL2) then return(false); fi;
  #
 
  vss2:=basisrec2.vss;
  vss2dual:=List(vss2, vs->vs*GramL2);
  getvs:=function(nrm)
    return(vss2[SinglePosition(basisrec2.nrmsset, nrm)]);
  end;
  getvsdual:=function(nrm)
    return(vss2dual[SinglePosition(basisrec2.nrmsset, nrm)]);
  end;
  #
  basis1:=basisrec1.basis;
  fingerprints1:=basisrec1.fingerprints;
  basis1inv:=InverseMat(basis1);
  newGram1:=basisrec1.newGram;
  #
  totalflag:=false;
  thetg:=[];
  #
  extend:=function(psol)
    #
    local leng, tnrm, tintnumbs, tvsdual, tvs, tvdual, pos,
    candidates, psoltvdual, tv;
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
    tnrm:=basisrec1.basisnrms[leng+1];
    if  leng=0 then 
      candidates:=getvs(tnrm);
    else 
      candidates:=[];
      tintnumbs:=List([1..leng], kk->newGram1[leng+1][kk]);
      tvsdual:=getvsdual(tnrm);
      tvs:=getvs(tnrm);
      pos:=0;
      for tvdual in tvsdual do 
        pos:=pos+1;
        psoltvdual:=psol*tvdual;
        if psoltvdual=tintnumbs then 
          Add(candidates, tvs[pos]);
        fi;
        #
        # this part cannot be changed to elif, because tintnumbs may be a zero vector
        #
        if  psoltvdual=-tintnumbs then 
          Add(candidates, -tvs[pos]);
        fi;
      od;
    fi;
    #
    if Length(candidates)<>fingerprints1[leng+1] then return(); fi;
    #
    for tv in candidates do 
      Add(psol, tv);
      extend(psol);
      if Remove(psol)<>tv then beep(776611); fi;
      if totalflag then return(); fi;
    od;
    #
  end;
  #
  extend([]);
  if totalflag then return(thetg);
  else return(fail);
  fi;
  #
end;



#####