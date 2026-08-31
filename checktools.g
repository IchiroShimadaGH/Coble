# Read("checktools.g");

Read("paperaffineConway/affineConwayCompdata.txt");

uumat:=[[0,1/2], [1/2,0]];
vvmat:=[[1,1/2], [1/2, 1]];
aamat:=[[1/2]];
bbmat:=[[3/2]];

Read("AutDiscfByStabilizerChain.g");
Read("OGLat.g");
Read("ImKerOLtoOqL.g");

readdata("isompairs");

isisomuvab:=function(uvab1, uvab2)
  local uvab1s, pos, tuvab, tp, tpp, dduvab, nnuvab, xx, 
  iseven1, iseven2, leng1, leng2;
  if uvab1=uvab2 then return(true); fi;
  #
  leng1:=[2,2,1,1]*uvab1;
  leng2:=[2,2,1,1]*uvab2;
  #
  if leng1<>leng2 then return(false); fi;
  #
  iseven1:=(uvab1[3]=0) and (uvab1[4]=0);
  iseven2:=(uvab2[3]=0) and (uvab2[4]=0);
  #
  if iseven1<>iseven2 then return(false); fi;
  #
  #
  uvab1s:=[uvab1];
  pos:=1;
  while pos<=Length(uvab1s) do 
    tuvab:=uvab1s[pos];
    for tp in  isompairs do
      tpp:=tp;
      dduvab:=tuvab-tpp[1];
      if Set(dduvab, xx->xx>=0)=[true] then 
        nnuvab:=dduvab+tpp[2];
        if nnuvab=uvab2 then return(true);fi;
        if not nnuvab in uvab1s then Add(uvab1s, nnuvab); fi;
      fi;
       tpp:=Reversed(tp);
       dduvab:=tuvab-tpp[1];
       if Set(dduvab, xx->xx>=0)=[true] then 
        nnuvab:=dduvab+tpp[2];
        if nnuvab=uvab2 then return(true);fi;
        if not nnuvab in uvab1s then Add(uvab1s, nnuvab); fi;
      fi;
    od;
    pos:=pos+1;
  od;
  return(false);
end;




uvabvect:=function(discrec)
  local ntqfrec, ntqf, xx, tuvab;
  if Set(discrec.discg)<>[2] then beep(52875872); fi;
  ntqfrec:=Normalform_fqf(discrec.discg, discrec.discf);
  ntqf:=ntqfrec.decompdata[1].pfqf;
  tuvab:=List([uumat, vvmat, aamat, bbmat], xx->OccurNumb(ntqf, xx));
  return(tuvab);
end;

uvab2qf:=function(tuvabvect)
  local dmatslist, uvabmatslist, ii, kk;
  dmatslist:=[];
  uvabmatslist:=[uumat, vvmat, aamat, bbmat];
  ii:=0;
  for kk in tuvabvect do
    ii:=ii+1;
    Append(dmatslist, List([1..kk], xx->uvabmatslist[ii]));
  od;
  return(DiagonalMats(dmatslist));
end;

#AutDiscfByStabilizerChain:=function(inibiis, bigmat)

uvabOqLOrder:=function(tuvabvect)
  local tqf, leng, vs, jj, inibiis, tv, bigmat, tsize;
  tqf:=uvab2qf(tuvabvect);
  leng:=[2,2,1,1]*tuvabvect;
  vs:=Cartesian(List([1..leng], jj->[0,1]));
  inibiis:=List(IdentityMat(leng), tv->SinglePosition(vs, tv));
  bigmat:=CopyNormalDiscf(TMTTmult(vs, tqf));
  tsize:=AutDiscfByStabilizerChain(inibiis, bigmat).order;
  return(tsize);
end; 


###########################3

OrbitdDecompByPerms:=function(nn, perms)
  #
  # perms is a list of permutations of [1..nn].
  # orbs is the list of orbits of the action of Group(permutations) on [1..nn].
  #
  local orbs, dones, tg, ii, jj, cc, torb, jjtg, forbs, beep;
  #
  beep:=function(beepnumb)
    localbeep("OrbitdDecompByPerms", beepnumb); Error();
  end;
  #
  orbs:=[];
  dones:=[];
  #
  for ii in [1..nn]  do 
    if ii in dones then continue; fi;
    torb:=[ii];
    cc:=0;
    for jj in torb do 
      cc:=cc+1;
      for tg in perms do 
        jjtg:=jj^tg;
        if not jjtg in torb then Add(torb, jjtg); fi;
      od;
    od;
    if cc<>Length(torb) then beep(91919); fi;
    Add(orbs, torb);
    Append(dones, torb);
  od;
  #
  Sort(dones);
  if dones<>[1..nn] then beep(71967981); fi;
  orbs:=List(orbs, Set);
  forbs:=Flat(orbs);
  Sort(forbs);
  if forbs<>[1..nn] then beep(71955981); fi;
  #
  return(orbs);
  #
end;




####

weyl0:=MakeVectei(26, 1);


toLeechRoot:=function(tlambda)
  local nlambda, tr;
  nlambda:=tlambda*GramLeech*tlambda;
  tr:=List(tlambda);
  Add(tr, 1, 1);
  Add(tr, nlambda/2-1, 1);
  if tr*GramL26*tr<>-2 then beep(51592952); fi;
  return(tr);
end;

toOGL26:=function(tgamma)
  local trs, tg, ttau, newtrs, tgg, w0, tvs, tv;
  tvs:=[];
  trs:=[];
  while Length(tvs)<40 or Rank(trs)<26  do 
    tv:=RandomVectFromL(24, [-2,-1,0,1,2]);
    Add(tvs, tv);
    Add(trs, toLeechRoot(tv));
  od;
  tg:=tgamma[1];
  ttau:=tgamma[2];
  newtrs:=List(tvs, tv->toLeechRoot(tv*tg+ttau));
  #AXisBwithFail:=function(A, B)
  tgg:=AXisBwithFail(trs, newtrs);
  if tgg=fail then beep(6199611); fi;
  if weyl0*tgg<>weyl0 then beep(581825); fi;
  if TMTTmult(tgg, GramL26)<>GramL26 then beep(919211); fi;
  return(tgg); 
end;


#########################


OqLtoOqLPerm:=function(tqg, discg, vs)
  local poss, vstqg, ttv, ii, leng;
  leng:=Length(discg);
  vstqg:=List(vs*tqg, ttv->List([1..leng], ii-> (ttv[ii] mod discg[ii])));
  poss:=List(vstqg, ttv->SinglePosition(vs, ttv));
  return(PermList(poss));
end;


#########################

MakeQF:=function(basisLeechm, tau24)
  #
  # QF(x) is xQx+2Lx+c.
  # The result is the quadratic form  x ->nn*((x-tau/2)^2-2) =
  # = nn (xQx-2*(tau/2)*Q*x +(tau/2)*Q*(tau/2)-2).
  # in the dual basis of Leechm.
  # htauindual is (tau/2) in the dual basis.
  # The integer nn is introduced to make everything in intger entries.
  #
  local tau, GramLeechm, htau, htauindual, Q,L,c,aas,nn, QFrec, GramLeechmdual, tn;
  #
  GramLeechm:=TMTTmult(basisLeechm, GramLeech);
  tn:=Length(GramLeechm);
  if Length(tau24)<>24 then beep(69698168); fi;
  tau:= SolutionMat(basisLeechm, tau24);
  if tau=fail then beep(661681); fi;
  if Length(tau)<>tn then beep(6692692); fi;
  if SignatureQ(GramLeechm)<>[tn, tn, 0] then beep(696991); fi;
  #
  GramLeechmdual:=InverseMat(GramLeechm);
  htau:=tau/2;
  htauindual:=htau*GramLeechm;
  Q:=GramLeechmdual;
  L:=-htauindual*Q;
  c:=htauindual*Q*htauindual-2;
  aas:=[];
  aas:=Concatenation(Flat(Q), aas);
  aas:=Concatenation(L, aas);
  Add(aas, c);
  nn:=Llcm(List(aas, DenominatorRat));
  Q:=nn*Q;
  L:=nn*L;
  c:=nn*c;
  if not IsIntMat(Q) then beep(194491); fi;
  if not IsIntVect(L) then beep(124491); fi;
  if not IsInt(c) then beep(144491); fi;
  QFrec:=rec(
    scalar:=nn,
    Q:=Q,
    L:=L,
    c:=c
  );
  return(QFrec);
end;


MakeXirec:=function(basisLeechm, ttau24)
  #
  local GramLeechm, GramLeechmdual, Xi, tilXi, ttask, tQ, tL, tc, QFrec, scalar,
  bdr, Xirec, ttau, xx;
  #
  GramLeechm:=TMTTmult(basisLeechm, GramLeech);
  GramLeechmdual:=InverseMat(GramLeechm);
  ttau:= SolutionMat(basisLeechm, ttau24);
  if ttau=fail then beep(64441); fi;
  #
  Xi:=[];
  tilXi:=[];
  bdr:=0;
  #
  QFrec:=MakeQF(basisLeechm, ttau24);
  scalar:=QFrec.scalar;
  tQ:=QFrec.Q;
  tL:=QFrec.L;
  tc:=QFrec.c;
  bdr:=0;
  #
  ttask:=function(sol)
    local tv, tx, ta, tb, xi, tilxi;
    tb:=sol*tQ*sol+2*tL*sol+tc;
    tv:=sol*GramLeechmdual;
    tx:=tv-ttau/2;
    ta:=tx*GramLeechm*tx;
    if (ta-2)*scalar<>tb then beep(5778); fi;
    if ta>2 then beep(89692); 
    elif ta<2 then 
      xi:=(List(tx)+ttau/2)*basisLeechm;
      tilxi:=tx*basisLeechm;
      Add(tilXi, tilxi);
      Add(Xi, xi);
    elif ta=2 then bdr:=bdr+1;
    else beep(66871); 
    fi;
    return(true);
  end;
  #
  #LLLIntVectsQF:=function(Q, L, c, isES, task)
  #
  LLLIntVectsQF(tQ, tL, tc, false, ttask);
  #
  if not IsEqualSet(tilXi, -tilXi) then beep(669629); fi;
  if Set(tilXi, xx->xx*GramLeech*xx<2)<>[true] then beep(229862); fi;
  if tilXi<>List(Xi, xx->xx-ttau24/2) then beep(339862); fi;
  #
  Xirec:=rec(
    tau:=ttau,
    Xi:=Xi, 
    tilXi:=tilXi,
    QF:=QFrec, 
    bdr:=bdr
  );
  return(Xirec);
end;



#########################

IsSigma12:=function(tGram)
  local ii, svsrec;
  #
  if not IsIntMat(tGram) then beep(686961); fi;
  if TransposedMat(tGram)<>tGram then beep(514276); fi;
  if SignatureQ(tGram)<>[12, 12, 0] then beep(716242); fi;
  if DeterminantIntMat(tGram)<>1 then beep(551121); fi;
  if Set([1..12], ii->(tGram[ii][ii] mod 2)=0)=[true] then beep(716221); fi;
  #
  svsrec:=ShortestVectors(tGram,  1);
  if svsrec.vectors<>[] then beep(9898981); fi;
  #
  return(true); 
end;

IsE8:=function(tGram)
  local ii, svsrec;
  #
  if not IsIntMat(tGram) then beep(9654441); fi;
  if TransposedMat(tGram)<>tGram then beep(569876); fi;
  if SignatureQ(tGram)<>[8, 8, 0] then beep(52254171); fi;
  if DeterminantIntMat(tGram)<>1 then beep(5511); fi;
  if Set([1..8], ii->(tGram[ii][ii] mod 2)=0)<>[true] then beep(712221); fi;
  #
  svsrec:=ShortestVectors(tGram,  2);
  if Length(svsrec.vectors)<>120 then beep(229981); fi;
  #
  return(true); 
end;


IsUnimodHyp10:=function(tGram, parity)
  local ii, diagps;
  #
  if not IsIntMat(tGram) then beep(6113386961); fi;
  if TransposedMat(tGram)<>tGram then beep(5122444276); fi;
  if SignatureQ(tGram)<>[10, 1, 9] then beep(7164455242); fi;
  if DeterminantIntMat(tGram)<>-1 then beep(555511251); fi;
  diagps:=Set([1..10], ii->(tGram[ii][ii] mod 2)=0);
  if parity=0 then 
    if diagps<>[true] then beep(71556221); fi;
  elif parity=1 then
    if diagps=[true] then beep(445256221); fi; 
  else beep(3456765434567);
  fi;
  #
  #
  return(true); 
end;

#########################

ChooseRandomWord:=function(n)
  local tw;
  if not n in [0,8,12,16,24] then beep(71711); fi;
  if n=0 then return([]);
  elif n=24 then return(theOmega);
  else 
    while true do
      tw:=Random(GolayWords);
      if Length(tw)=n then return(tw); fi;
    od;
  fi;
  return();
end;

LeechBasisInv:=InverseMat(LeechBasis);
theid24:=IdentityMat(24);
theid26:=IdentityMat(26);

MakeVep:=function(word)
  local tgb, poss, aa, tpos, tg;
  tgb:=List(theid24);
  poss:=List(word, aa->SinglePosition(theOmega, aa));
  for tpos in poss do tgb[tpos]:=-tgb[tpos]; od;
  tg:=LeechBasis*tgb*LeechBasisInv;
  if not IsIntMat(tg)then beep(761221); fi;
  if TMTTmult(tg, GramLeech)<>GramLeech then beep(55112); fi;
  return(tg);
end;

ChooseRandomTau:=function(vep)
  local basism, rv, tau;
  basism:=NullspaceIntMat(theid24+vep);
  rv:=RandomVectFromL(Length(basism), [-2,-1,0,1,2]);
  tau:=rv*basism;
  if tau*vep<>-tau then beep(611761); fi;
  return(tau);
end;

TypeTau:=function(vep, tau)
  local impi, tv, aa;
  impi:=theid24+vep;
  tv:=SolutionIntMat(impi, tau);
  if tv<>fail then return("zero"); fi;
  aa:=(1/2)*(tau*GramLeech*tau);
  aa:=aa mod 2;
  if aa=1 then return("120"); 
  else return("135"); 
  fi;
end;

ChooseRandomTauType:=function(nn, vep, type)
  local basism, impi, cc, rv, tv, aa, tau;
  basism:=NullspaceIntMat(theid24+vep);
  impi:=theid24-vep;
  if not type in ["zero", "120", "135"] then beep(6891691); fi;
  if nn<>16 and type<>"zero" then beep(919192); fi;
  cc:=0;
  if type="zero" then 
    rv:=RandomVectFromL(Length(impi), [-2,-1,0,1,2]);
    tau:=rv*impi;
    return(tau);
  fi;
  while true do 
    rv:=RandomVectFromL(Length(basism), [-2,-1,0,1,2]);
    tau:=rv*basism;
    if tau*vep<>-tau then beep(611761); fi;
    tv:=SolutionIntMat(impi, tau);
    if tv=fail then 
      aa:=(1/2)*(tau*GramLeech*tau);
      aa:=aa mod 2;
      if aa=1 then if type="120"  then return(tau); fi; 
      elif aa=0 then if type="135"  then return(tau); fi; 
      fi;
    fi;
    cc:=cc+1;
    if cc=10000 then beep(182871); fi;
  od;
end;


NewChooseRandomTauType:=function(nn, vep, type)
  local basism, impi, cc, rv, tv, aa, tau;
  basism:=NullspaceIntMat(theid24+vep);
  impi:=theid24-vep;
  if not type in ["zero", "120", "135"] then beep(6891691); fi;
  if nn<>16 and type<>"zero" then beep(919192); fi;
  cc:=0;
  if type="zero" then 
    rv:=RandomVectFromL(Length(impi), [-2,-1,0,1,2]);
    tau:=rv*impi;
    return(tau);
  fi;
  while true do 
    rv:=RandomVectFromL(Length(basism), [-2,-1,0,1,2]);
    tau:=rv*basism;
    if tau*vep<>-tau then beep(611761); fi;
    tv:=SolutionIntMat(impi, tau);
    if tv<>fail then 
      if type="zero" then return(tau); fi;
    else 
      aa:=(1/2)*(tau*GramLeech*tau);
      aa:=aa mod 2;
      if aa=1 then if type="120"  then return(tau); fi; 
      elif aa=0 then if type="135"  then return(tau); fi; 
      fi;
    fi;
    cc:=cc+1;
    if cc=100000000 then beep(182871); fi;
  od;
end;



weyl0:=MakeVectei(26, 1);


toLeechRoot:=function(tlambda)
  local nlambda, tr;
  nlambda:=tlambda*GramLeech*tlambda;
  tr:=List(tlambda);
  Add(tr, 1, 1);
  Add(tr, nlambda/2-1, 1);
  if tr*GramL26*tr<>-2 then beep(51592952); fi;
  return(tr);
end;

GetGamma:=function(vep, tau)
  local trs, tg, ttau, newtrs, tgg, w0, tvs, tv;
  tvs:=[];
  trs:=[];
  while Length(tvs)<40 or Rank(trs)<26  do 
    tv:=RandomVectFromL(24, [-3, -2,-1,0,1,2, 3]);
    Add(tvs, tv);
    Add(trs, toLeechRoot(tv));
  od;
  newtrs:=List(tvs, tv->toLeechRoot(tv*vep+tau));
  #AXisBwithFail:=function(A, B)
  tgg:=AXisBwithFail(trs, newtrs);
  if tgg=fail then beep(6199611); fi;
  if not IsIntMat(tgg) then beep(312131); fi;
  if weyl0*tgg<>weyl0 then beep(581825); fi;
  if TMTTmult(tgg, GramL26)<>GramL26 then beep(919211); fi;
  return(tgg); 
end;


GetGBrec:=function(vep, tau)
  local basisp, basism, GBrec, gamma;
  gamma:=GetGamma(vep, tau);
  basisp:=NullspaceIntMat(theid26-gamma);
  basism:=NullspaceIntMat(theid26+gamma);
  if CokerTorsion(basisp)<>[] then beep(6252414); fi;
  if CokerTorsion(basism)<>[] then beep(5452414); fi;
  GBrec:=rec(
    vep:=vep,
    tau:=tau, 
    gamma:=gamma,
    basisp:=basisp,
    basism:=basism,
    Gramp:=TMTTmult(basisp, GramL26),
    Gramm:=TMTTmult(basism, GramL26)
  );
  return(GBrec);
end;

theg8:=174182400;
theg12:=2^11*Factorial(12);


CheckGBrec:=function(n, type, GBrec)
  #
  local beep, Gramp, Gramm, svsrec, grec, erec, nf,
  mOGrec, discv, discvv, imqsize, Kerqrec,  basisLeechm,
  tilXi, tilXXi, col1, col2, xx, vv;
  #
  beep:=function(beepnumb)
    localbeep("CheckGBrec", beepnumb); Error();
  end;
  #
  Gramp:=GBrec.Gramp;
  Gramm:=GBrec.Gramm;
  if SignatureQ(Gramp)<>[26-n, 1, 25-n] then beep(3881211); fi;
  if SignatureQ(Gramm)<>[n, 0, n]  then beep(1881211); fi;
  #
  svsrec:=ShortestVectors(-Gramm, 4);
  if Set(svsrec.norms)<>[4] then beep(919911); fi;
  nf:=2*Length(svsrec.vectors);
  #
  mOGrec:=OGLat(-Gramm);
  if NewGensSize(-Gramm, mOGrec.totalgens)<>mOGrec.order then beep(88181); fi;
  #
  Kerqrec:=KerOLtoOqL(-Gramm);
  imqsize:=ImOLtoOqL(mOGrec.totalgens, -Gramm);
  CheckKerq(-Gramm, Kerqrec.gens);
  if NewGensSize(-Gramm, Kerqrec.gens)<>Kerqrec.size then beep(66181); fi;
  if Kerqrec.size*imqsize<>mOGrec.order then beep(318181); fi;
  #
  basisLeechm:=NullspaceIntMat(theid24+GBrec.vep);
  tilXi:=MakeXirec( basisLeechm, GBrec.tau).tilXi;
  #
  discv:=uvabvect(DiscriminantForm(Gramm));
  #
  if n=8 then 
    if type<>"zero" then beep(888121); fi;
    erec:=erecs[1];
    grec:=gammarecs[1];
    if nf<>240 then beep(91912); fi;
    discvv:=[4,0,0,0];
    if mOGrec.order<> 4*theg8 then beep(99991); fi;
    if not IsE8(-Gramm/2) then beep(818181); fi;
  elif n=12 then 
    if type<>"zero" then beep(338121); fi;
    erec:=erecs[2];
    grec:=gammarecs[2];
    if nf<>264 then beep(93312); fi;
    discvv:=[0,0,12,0];
    if mOGrec.order<> theg12 then beep(99991); fi;
    if not IsSigma12(-Gramm/2) then beep(818181); fi;
  elif n=16 then 
     erec:=erecs[3];
    if type="zero" then 
      grec:=gammarecs[3];
      if nf<>4320 then beep(73312); fi;
      discvv:=[4,0,0,0];
      if mOGrec.order<> 512*theg8 then beep(99991); fi;
    elif type="120" then 
      grec:=gammarecs[4];
      if nf<>2016 then beep(53312); fi;
      discvv:=[0,0,5,5];
      if mOGrec.order<> 743178240 then beep(9229991); fi;
      if not IsUnimodHyp10(Gramp/2, 1) then beep(77162); fi;
    elif type="135" then 
      grec:=gammarecs[5];
      if nf<>2272 then beep(13312); fi;
       discvv:=[5,0,0,0];
      if mOGrec.order<> 21139292160 then beep(955991); fi;
      if not IsUnimodHyp10(Gramp/2, 0) then beep(755262); fi;
    else 
      beep(998121); 
    fi;
  else beep(621189);
  fi;
  #
  if not isisomuvab(discv, discvv) then beep(576457); fi;
  if mOGrec.order<>grec.minusrec.OGrec.size then beep(93911); fi;
  if Kerqrec.size<>grec.minusrec.Kerqrec.size then beep(79911); fi;
  #
  col1:=Collected(List(tilXi, xx->xx*GramLeech*xx));
  tilXXi:=List(grec.plusrec.Xitau, vv->vv-grec.tau/2);
  col2:=Collected(List(tilXXi, xx->xx*GramLeech*xx));
  if col1<>col2 then beep(162626); fi;
  #
  return(true);
end;

#################