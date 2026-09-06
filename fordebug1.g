NewOGLat:=function(GramL, basisrec)
  local beep,  nn, LLLrec, basis, newGram, basisnrmslist,
  ii, maxnrm, svrec, nrmsset, pvposss, pvlists, aa, poss, pos, getpvlist,
  tbasis, tbasisdual, fingerprints, lllist, minll, inipos, 
  initb, OGrec, getfp,  tfps, tbs, tb, minfp, newtb,
  basisinv, basisdual, basisnrms, reachflag, thetg, simpleextend,
  getorb, candidates1, getpermcan1, stabrecs, stabrec, doneposs,
  gengs, genperms, trueposs, neworb, tv, candidates, getpermcan, tperm,
  inipsol, kk, GramLinv, pvduallists, getpvduallist, tpvlist,
  getorbsunion, falseseeds, trueseeds, falseposs, thelevel, inipsoldual,
  tintnumbs, tpv, jj, intnumbss, totalgens, lengcan, lengcan1, tblist, fblist,
  getorbblist, getorbsunionblist;
  #
  beep:=function(beepnumb)
    localbeep("OGLat", beepnumb); Error();
  end;
  #
  nn:=Length(GramL);
  basis:=basisrec.basis;
  newGram:=TMTTmult(basis, GramL);
  basisnrmslist:=List([1..nn], ii->newGram[ii][ii]);
  maxnrm:=basisrec.maxnrm;
  nrmsset:=basisrec.nrmsset;
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
  pvlists:=List(basisrec.vss, vs->List(vs,  inipos));
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
  # We permute basis to make fingerprints smaller.
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
  od;
  #
  # The end of "We permute basis to make fingerprints smaller".
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
  #
  # level 1 is done
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
    size:=2*Product(List(stabrecs, stabrec->stabrec.order)),
    totalgens:=totalgens
  );
  #
  return(OGrec);
  #
end;

IsIsomBasisRecs:=function(basisrec1, basisrec2)
  if basisrec1.
end;



####