#Read("OGLat.g");

#PleskenSouvignier

OGLat:=function(arg)
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
    localbeep("OGLat", beepnumb); Error();
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
# OGrecL1 is the existing, unmodified result of OGLat(GramL1).
# Positive definite integral Gram matrices of positive rank are assumed.
# Returns fail iff the lattices are not isomorphic; otherwise returns T with
#   T*GramL2*TransposedMat(T) = OGrecL1.Gram and Det(T) = +/-1.
# In row coordinates the isomorphism L1 -> L2 is v -> v*T.
# Neither OGLat nor an automorphism computation for L2 is called.

IsomLats:=function(OGrecL1, GramL2)
  #
  local GramL1, n, basis1, basis1inv, sourceGram, norms,
  red2, basis2, targetGram, sv, shells, duals, fingerprints,
  positive, search, result, first, getshell, beep, ss, vv;
  #
  beep:=function(beepnumb)
    localbeep("IsomLats", beepnumb); Error();
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
  # These lengths recover the fingerprints from the existing OGLat record.
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