#Read("newtools.g");
#


readdata("GolayVects");
readdata("GolayWords");
readdata("theOmega");

readdata("GramU");
readdata("GramE8");
readdata("GramLeech");
readdata("GramL26");
readdata("weyl0");
readdata("LeechBasis");

readdata("SizeOfG");


AutDiscfByStabilizerChain:=function(inibiis, bigmat)
  #
  local beep, nn, targetmat, ii, jj, orbits, endmark, task, autrec, examples,
  yy, kk, tv;
  #
  # We have made the list of all the elements of discg.
  # We have also made bigmat from  the list.
  # inibiis is the basis of discg given by the positions in the list.
  #
  beep:=function(beepnumb)
    localbeep("AutDiscfByStabilizerChain", beepnumb);Error();
  end;
  #
  nn:=Length(inibiis);
  if bigmat<>TransposedMat(bigmat) then beep(68215873); fi;
  targetmat:=List(inibiis, ii->List(inibiis, jj->bigmat[ii][jj]));
  orbits:=List([1..nn], ii->[]);
  examples:=List([1..nn], ii->[]);
  endmark:="end";
  #
  task:=function(piis, leng)
    #
    local tmm, xx,mm,targetv, kk, newpii, newpiis, rmm, pkk;
    #
    if leng=nn then
      tmm:=0;
      for xx in piis do
        if xx<>inibiis[tmm+1] then
          break; #from for xx in piis do
        else
          tmm:=tmm+1;
        fi;
      od;
      mm:=tmm;
      for kk in [1..mm+1] do
        if kk<=nn then
          pkk:=piis[kk];
          if not pkk in orbits[kk] then
            Add(orbits[kk], pkk);
            Add(examples[kk], List(piis));
          fi;
        elif kk=nn+1 then
          if piis<>inibiis then beep(522457); fi;
        else
          beep(6761871);
        fi;
      od;
      #
      return(mm);
      #
    else #if leng=nn then
      #
      targetv:=List([1..leng+1], kk->targetmat[leng+1][kk]);
      newpii:=0;
      for tv in bigmat do
        newpii:=newpii+1;
        newpiis:=CopyAdd(piis, newpii);
        if targetv=List(newpiis, kk->tv[kk]) then
          rmm:=task(newpiis, leng+1);
          if rmm<>endmark and rmm<leng then
            return(rmm);
          fi;
        fi;
      od;
    fi;
    #
    return(endmark);
  end;#if leng=nn then
  #
  task([], 0);
  #
  for kk in [1..nn] do
    if not IsDuplicateFree(orbits[kk]) then beep(698621); fi;
    if orbits[kk]<>List(examples[kk], yy->yy[kk]) then
      beep(8758751);
    fi;
  od;
  #
  autrec:=rec(
    orbits:=orbits,
    examples:=examples,
    order:=Product(List(orbits, Length))
  );
  #
  return(autrec);
end;


#######

New2FindAnIsomDiscf:=function(discgf1, discgf2)
  local beep, discf1, discg1, discf2, discg2,  big2,
  iterg2, ii, size, vects2, xx, xxdual, yy, aa,
  bone2, tisom, leng2, fflag, findanisom, jj, leng1, ff, jj2flag;
  #
  beep:=function(beepnumb)
    localbeep("New2FindAnIsomDiscf", beepnumb); Error();
  end;
  #
  discg1:=discgf1[1];
  discg2:=discgf2[1];
  size:=Product(discg1);
  if size<>Product(discg2) then return(fail); fi;
  #
  discf1:=discgf1[2];
  discf2:=discgf2[2];
  big2:=NullMat(size, size);
  iterg2:=IteratorOfCartesianProduct(List(discg2, ff->[0..ff-1]));
  vects2:=[];
  ii:=0;
  for xx in iterg2 do
    ii:=ii+1;
    xxdual:=xx*discf2;
    big2[ii][ii]:=modtZ(xxdual*xx);
    jj:=0;
    for yy in vects2 do
      jj:=jj+1;
      aa:=modZ(yy*xxdual);
      big2[ii][jj]:=aa;
      big2[jj][ii]:=aa;
    od;
    Add(vects2, xx);
  od;
  #
  bone2:=DiagonalMat(discg2);
  tisom:=fail;
  leng1:=Length(discg1);
  fflag:=false;
  #
  #Printn(size);
  #
  findanisom:=function(tempiis)
    #
    local ttisom, tbone, newii, tv1, newaa, jj2, kk, tii, big2jj2, newtempiis;
    #
    #
    if Length(tempiis)=leng1 then
      ttisom:=List(tempiis, kk->vects2[kk]);;
      tbone:=StackMats(bone2, ttisom);
      if CokerTorsion(tbone)=[] then
        tisom:=List(ttisom);
        fflag:=true;
      fi;
    else
      newii:=Length(tempiis)+1;
      tv1:=discf1[newii];
      newaa:=tv1[newii];
      for jj2 in [1..size] do
        if big2[jj2][jj2]=newaa then
          jj2flag:=true;
          big2jj2:=big2[jj2];
          kk:=0;
          for tii in tempiis do
            kk:=kk+1;
            if big2jj2[tii]<>tv1[kk] then
              jj2flag:=false;
              break; #from for tii in tempiis do
            fi;
          od;
          if jj2flag then
            newtempiis:=CopyAdd(tempiis, jj2);
            findanisom(newtempiis);
            if fflag then
              break;
            fi;
          fi;
          if fflag then
            break;#from for jj2 in [1..size] do
          fi;
        fi;
      od;
    fi;
  end;
  #
  findanisom([]);
  #
  if fflag then
    if discf1<>CopyNormalDiscf(TMTTmult(tisom, discf2)) then
      beep(667267);
    fi;
  else
    if tisom<>fail then beep(565167); fi;
  fi;
  #
  return(tisom);
end;

NewRandomEvenLatFromL:=function(dim, L)
	local cc, ii, jj, aa, G;
	cc:=0;
	while true do
		G:=NullMat(dim, dim);
    for ii in [1..dim] do
      G[ii][ii]:=2*Random(L);
      for jj in [1..ii-1] do
        aa:=Random(L);
        G[ii][jj]:=aa;
        G[jj][ii]:=aa;
      od;
    od;
		if DeterminantIntMat(G)<>0 then return(G); fi;
		cc:=cc+1;
		if cc>10 then return(fail); fi;
	od;
	return(fail);
end;


NewRandomEvenLatFromLBdd:=function(dim, L, bdd)
	local cc, ii, jj, aa, G, absdet;
	cc:=0;
	while true do
		G:=NullMat(dim, dim);
    for ii in [1..dim] do
      G[ii][ii]:=2*Random(L);
      for jj in [1..ii-1] do
        aa:=Random(L);
        G[ii][jj]:=aa;
        G[jj][ii]:=aa;
      od;
    od;
    absdet:=AbsInt(DeterminantIntMat(G));
		if absdet>0 and absdet<=bdd then
      return(G);
    fi;
		cc:=cc+1;
		if cc>10 then return(fail); fi;
	od;
	return(fail);
end;


#############

#
MakeOGrec:=function(Gram)
  local beep, OGrec, torder, tRs, tbs, kk, tR, tg,
  bkk, bkktR, partbs, inirec, vg, n;
  #
  beep:=function(beepnumb)
    localbeep("MakeOGrec", beepnumb); Error();
  end;
  #
  n:=Length(Gram);
  if SignatureQ(Gram)<>[n,n,0] then beep(5476247634); fi;
  inirec:=OGPosLat(Gram, 0, 1);
  tRs:=List(inirec.orbs, torb->List(torb, vg->vg[2]));
  tbs:=List(inirec.basis);
  #
  #
  torder:=2*Product(List(tRs, Length));
  #
  # check
  #
  if torder<>inirec.order then beep(7652765); fi;
  if AbsInt(DeterminantIntMat(tbs))<>1 then beep(692521); fi;
  #
  kk:=0;
  for tR in tRs do
    kk:=kk+1;
    if Set(tR, tg->TMTTmult(tg, Gram)=Gram)<>[true] then
      beep(669769198);
    fi;
    bkk:=tbs[kk];
    bkktR:=List(tR, tg->bkk*tg);
    if not IsDuplicateFree(bkktR) then
      beep(663333198);
    fi;
    if kk=1 then
      if Intersection(bkktR, -bkktR)<>[] then beep(4465456); fi;
    fi;
    if kk>1 then
      partbs:=List([1..kk-1], jj->tbs[jj]);
      if Set(tR, tg->partbs*tg=partbs)<>[true] then
        beep(663369198);
      fi;
    fi;
  od;
  #
  #
  OGrec:=rec(
    Gram:=List(Gram),
    bs:=tbs,
    Rs:=tRs,
    order:=torder
  );
  #
  return(OGrec);
end;


WriteInOGrec:=function(gs, OGrec)
  local beep, n, tg, writeInOGrec, Rs, orbs, Idn, bs, ttg, kk;
  #
  beep:=function(beepnumb)
    localbeep("WriteInOGrec", beepnumb); Error();
  end;
  #
  n:=Length(OGrec.Gram);
  Idn:=IdentityMat(n);
  bs:=OGrec.bs;
  Rs:=OGrec.Rs;
  orbs:=[];
  for kk in [1..n] do
    Add(orbs, List(Rs[kk], ttg->bs[kk]*ttg));
  od;
  #
  writeInOGrec:=function(tg)
    #
    local ttg, rs, b1ttg, pos, r1, rkk, kk, bkkttg;
    #
    ttg:=List(tg);
    rs:=[];
    #
    b1ttg:=bs[1]*ttg;
    if b1ttg in orbs[1] then
      pos:=SinglePosition(orbs[1], b1ttg);
      r1:=Rs[1][pos];
      Add(rs, r1);
      ttg:=ttg*InverseMat(r1);
    elif -b1ttg in orbs[1] then
      pos:=SinglePosition(orbs[1], -b1ttg);
      r1:=Rs[1][pos];
      Add(rs, -Idn);
      Add(rs, r1);
      ttg:=ttg*InverseMat(-r1);
    else beep(365415452);
    fi;
    #
    for kk in [2..n] do
      bkkttg:=bs[kk]*ttg;
      pos:=SinglePosition(orbs[kk], bkkttg);
      rkk:=Rs[kk][pos];
      Add(rs, rkk);
      ttg:=ttg*InverseMat(rkk);
    od;
    #
    if ttg<>Idn then beep(5875782); fi;
    #
    rs:=Reversed(rs);
    if Product(rs)<>tg then beep(435678); fi;
    #
    return(rs);
  end;
  #
  for tg in gs do
    writeInOGrec(tg);
  od;
  #
  return(true);
end;

#########

RandomMatFromLSmallDet:=function( r,  L, trial )
  local tM, tdet, tt, ttM, ttdet;
  tM:=fail;
  tdet:=infinity;
  for tt in [1..trial] do
    ttM:=RandomMatFromL(r, r, L);
    ttdet:=AbsInt(DeterminantIntMat(ttM));
    if ttdet=0 then continue; fi; #in for tt in [1..trial] do
    if ttdet<tdet then
      tM:=ttM;
      tdet:=ttdet;
      if tdet=1 then break; fi; #from for tt in [1..trial] do
    fi;
  od;
  return(tM);
end;

######


Makeepsilonrec:=function(inicword)
  local beep, erec, LeechBasisInv, cword, nn,
  Id24, cvect,ii,tZmat, aa, tmat,
  plusbasis, minusbasis, plusGram, minusGram,
  pig, plusOGrec, minusOGrec;
  #
  beep:=function(beepnumb)
    localbeep("Makeepsilonrec", beepnumb); Error();
  end;
  #
  cword:=Set(inicword);
  if not cword in GolayWords then beep(8558123); fi;
  LeechBasisInv:=InverseMat(LeechBasis);
  #
  nn:=Length(inicword);
  #
  if not nn in [8,12,16,24] then beep(999123); fi;
  Id24:=List(Idmat24);
  cvect:=List([1..24], kk->0);
  ii:=0;
  tZmat:=[];
  for aa in theOmega do
    ii:=ii+1;
    aa:=theOmega[ii];
    if aa in cword then
      Add(tZmat, -Id24[ii]);
      cvect[ii]:=1;
    else
      Add(tZmat, Id24[ii]);
    fi;
  od;
  #
  if not cvect in GolayVects then beep(8668123); fi;
  #
  tmat:=LeechBasis*tZmat*LeechBasisInv;
  if not IsIntMat(tmat) then beep(771234); fi;
  if TMTTmult(tmat, GramLeech)<>GramLeech then beep(88833); fi;
  if Order(tmat)<>2 then beep(881231); fi;
  #
  Id24:=List(Idmat24);
  plusbasis:=NullspaceIntMat(tmat-Id24);
  if Length(plusbasis)<>24-nn then beep(88811); fi;
  if nn<>24 then
    plusbasis:=HermiteNormalFormIntegerMat(plusbasis);
    plusGram:=TMTTmult(plusbasis, GramLeech);
  else
    plusbasis:=[];
    plusGram:=[];
  fi;
  #
  minusbasis:=NullspaceIntMat(tmat+Id24);
  if Length(minusbasis)<>nn then beep(82211); fi;
  minusbasis:=HermiteNormalFormIntegerMat(minusbasis);
  minusGram:=TMTTmult(minusbasis, GramLeech);
  if nn=24 then
    if minusbasis<>Id24 then beep(68827); fi;
  fi;
  #
  pig:=MisXB(Id24-tmat, minusbasis);
  if not IsIntMat(pig) then beep(78698); fi;
  if Rank(pig)<>nn then beep(99991); fi;
  Printn("Coker pig", CokerTorsion(pig));
  #
  if nn<>24 then
    Printn("plusOG  start");
    plusOGrec:=MakeOGrec(plusGram);
    Printn("plusOG done", plusOGrec.order);
    Printn("minusOG  start");
    minusOGrec:=MakeOGrec(minusGram);
    Printn("minusOG done", minusOGrec.order);
  else
    plusOGrec:="NULL";
    minusOGrec:="not available";
  fi;
  #
  erec:=rec(
    n:=nn,
    cword:=List(cword),
    cvect:=cvect,
    mat:=tmat,
    plusbasis:=plusbasis,
    minusbasis:=minusbasis,
    plusOGrec:=plusOGrec,
    minusOGrec:=minusOGrec,
    piepsilon:=pig
  );
  return(erec);
end;


#################


OqLtoOqLPerm:=function(tqg, discg, vs)
  local poss, vstqg, ttv, ii, leng;
  leng:=Length(discg);
  vstqg:=List(vs*tqg, ttv->List([1..leng], ii-> (ttv[ii] mod discg[ii])));
  poss:=List(vstqg, ttv->SinglePosition(vs, ttv));
  return(PermList(poss));
end;

uumat:=[[0,1/2], [1/2,0]];
vvmat:=[[1,1/2], [1/2, 1]];
aamat:=[[1/2]];
bbmat:=[[3/2]];

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


#############
# 2026/04/26
#############

latexFactorsInt:=function(nn)
	local absnn, str, ii, tp, concat, pnus, ps, xx;
	str:="";
	if nn=0 then str:="0"; return(str); fi;
	if nn<0 then str:="-"; fi;
	absnn:=AbsInt(nn);
  if absnn=1 then 
    Append(str, "1");
    return(str);
  fi;
  pnus:=Collected(FactorsInt(absnn));
  ps:=List(pnus, xx->xx[1]);
  if not IsSortedList(ps) then 
    beep(892982829);
  fi;
	ii:=0;
	for tp in pnus do
		if ii>0 then
			Append(str, " \cdot ");
		fi;
		concat:=Concatenation("{", String(tp[1]), "}");
		Append(str, concat);
		if tp[2]>1 then
			concat:=Concatenation("^{", String(tp[2]), "}");
			Append(str, concat);
		fi;
		ii:=ii+1;
	od;
	return(str);
end;



GramEn:=function(n)
  local tG, edges, ii, te;
  if n<4 then beep(34657687); fi;
  tG:=(-2)*IdentityMat(n);
  edges:=[[1,4]];
  for ii in [2..n-1] do
    Add(edges, [ii, ii+1]);
  od;
  for te in edges do
    tG[te[1]][te[2]]:=1;
    tG[te[2]][te[1]]:=1;
  od;
  return(tG);
end;

#######

#############
# 2026/05/10
#############

## from Task20260423a.g

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


############# Added om 2026/06/12 based on Task20260416a

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

#################

nowsize:=490497638400;
targetsize:=51231497335603200;

readdata("Sigma12rec");
readdata("OqGL12pluspermsIni");
readdata("GramLp14");
readdata("discrecLp14");
readdata("discvsLp14");

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


GetFramSigma12:=function(tGram)
  local svsrec, twovects, cc1, cc2, cc3, tv1, tv2, tv3,
  tv1dual, tv2dual, tv3dual, triangles, counter, xx, tfvs,
  frame, tv, tsum, dones, triangle, tpair, tvdual;
  #
  if not IsSigma12(tGram) then beep(55442); fi;
  svsrec:=ShortestVectors(tGram,  2);
  twovects:=Union(svsrec.vectors, -svsrec.vectors);
  triangles:=[];
  cc1:=0;
  for tv1 in twovects do
    cc1:=cc1+1;
    tv1dual:=tv1*tGram;
    cc2:=0;
    for tv2 in twovects do
      cc2:=cc2+1;
      if cc1>=cc2 or  tv2*tv1dual<>1 then continue; fi;# in for tv2 in twovects do
      tv2dual:=tv2*tGram;
      cc3:=0;
      for tv3 in twovects do
        cc3:=cc3+1;
        if cc2>=cc3 or tv3*tv1dual<>1 or tv3*tv2dual<>1 then continue; fi;# in for tv3 in twovects 
        tv3dual:=tv3*tGram;
        #
        counter:=0;
        #
        for xx in twovects do 
          if xx*tv1dual=1 and xx*tv2dual=1 and xx*tv3dual=1 then 
            counter:=counter+1;
          fi;
        od;
        #
        if counter=0 then 
          Add(triangles, [tv1, tv2, tv3]);
        else 
          if counter<>16 then beep(226981); fi;
        fi;
        #
      od;#for tv3 in twovects do
    od; #for tv2 in twovects do
  od;#for tv1 in twovects d
  #
  if Length(triangles)<>1760 then beep(413211); fi;
  #
  tfvs:=[];
  for triangle in triangles do
    tsum:=Sum(triangle);
    for tv in triangle do 
      Add(tfvs, tsum-2*tv);
    od; 
  od;
  #
  tfvs:=Set(tfvs);
  #
  if Length(tfvs)<>24 then beep(698691); fi;
  if not IsEqualSet(tfvs, -tfvs) then beep(881631); fi;
  #
  frame:=[];
  dones:=[];
  #
  for tv in tfvs do 
    if not tv in dones then 
      tvdual:=tv*tGram;
      if frame<>[ ] then 
        if Set(frame, tpair->tpair*tvdual )<>[[0,0]] then 
          beep(717621); 
        fi;
      fi;
      Add(frame, [tv, -tv]);
      Append(dones, [tv, -tv]);
    fi;
  od;
  #
  if Length(frame)<>12 then beep(338691); fi;
  if Set(frame, Length)<>[2] then beep(338331); fi;
  #
  return(frame);
  #
end;

FindIsomSigma12s:=function(Sigmarec, tGramSigma) 
  local Gram0, framebasis0inv, frame1, framebasis1, tg, tpos; 
  Gram0:=Sigmarec.Gram;
  framebasis0inv:=Sigmarec.framebasisinv;
  frame1:=GetFramSigma12(tGramSigma);
  framebasis1:=List(frame1, xx->xx[Random([1,2])]);
  tg:=framebasis0inv*framebasis1;
  if not (IsIntMat(tg) and TMTTmult(tg, tGramSigma)=Gram0) then 
    tpos:=Random([1..12]);
    framebasis1[tpos]:=-framebasis1[tpos];
    tg:=framebasis0inv*framebasis1;
    if not (IsIntMat(tg) and TMTTmult(tg, tGramSigma)=Gram0) then
      beep(443121); 
    fi; 
  fi;
  return(tg);
end;

GetnewsizeSigma12:=function(tUB)
  local newtGram,  onevects, ttg, tgL, tperm, ttgLperms, newttsize, tGG, sumtUB,
  orthorec;
  #
  #
  if TMTTmult(tUB, GramLp14)<>[[0,1],[1,0]] then beep(531212); fi;
  sumtUB:=Sum(tUB);
  if sumtUB[1]+sumtUB[2]<=0 then beep(98181); fi;
  orthorec:=OrthogonalCompRec(GramLp14, tUB);
  newtGram:=(-1/2)*orthorec.Gram;
  if not IsIntMat(newtGram) then beep(33969999); fi;
  if not IsSymmMat(newtGram) then beep(33922999); fi;
  if SignatureQ(newtGram)<>[12,12,0] then beep(53352999); fi;
  if DeterminantIntMat(newtGram)<>1 then beep(37165211); fi;
  if IsEvenLattice(newtGram) then beep(3666191); fi;
  onevects:=ShortestVectors(newtGram, 1).vectors;
  #
  if onevects=[] then 
    #
    ttg:=FindIsomSigma12s(Sigma12rec, newtGram);
    tgL:=Concatenation(tUB, ttg*orthorec.basis);
    if not IsIntMat(tgL) then beep(3251427); fi;
    if TMTTmult(tgL, GramLp14)<>GramLp14 then beep(7712221);fi;
    tperm:=OqLtoOqLPerm(OLtoOqL(tgL, discrecLp14), discrecLp14.discg, discvsLp14);
    ttgLperms:=List(OqGL12pluspermsIni);
    Add(ttgLperms, tperm);
    newttsize:=Size(Group(ttgLperms));
    if newttsize>nowsize then 
      if newttsize<>targetsize then beep(719981); fi;
      return("Sigma12, new");
    elif newttsize=nowsize then 
      return("Sigma12, not new");
    else beep(797981); 
    fi;
  else 
    if CokerTorsion(onevects)<>[] then beep(51551); fi;
    if Length(onevects)=4 then 
      if TMTTmult(onevects, newtGram)<>Idmat4 then beep(876871); fi;
      tGG:=OrthogonalCompRec(newtGram, onevects).Gram;
      if DeterminantIntMat(tGG)<>1 then beep(375211); fi;
      if not IsEvenLattice(tGG) then beep(366191); fi; 
      if ShortestVectors(tGG, 1).vectors<>[] then beep(995221); fi;
      if Length(ShortestVectors(tGG, 2).vectors)<>120  then beep(775221); fi;
      return("E8+I4");
    elif Length(onevects)=12 then 
      if TMTTmult(onevects, newtGram)<>Idmat12 then beep(876871); fi;
      return("I12");
    else beep(5875871); 
    fi;
  fi;
  #
  beep(56121);
  return();
  #
end; 

#
# Made and checked in 
# Read("Task20260622a.g");
#

readdata("LeechBasisinv");
readdata("M24gens");


M24toOGLeech:=function(tgperm)
  local tI, ttgmat, ii,  tgmat;
  tI:=Idmat24;
  ttgmat:=List([1..24], ii->tI[ii^tgperm]);
  tgmat:=LeechBasis*ttgmat*LeechBasisinv;
  if not IsIntMat(tgmat) then beep(77112); fi;
  if TMTTmult(tgmat, GramLeech)<>GramLeech then beep(11241); fi;
  return(tgmat);
end;

GolaytoOGLeech:=function(tword)
  local tI, ttgmat, poss, pos,  tgmat;
  poss:=Set(tword, xx->Position(theOmega, xx));
  ttgmat:=List(Idmat24);
  for pos in poss do
    ttgmat[pos]:=-ttgmat[pos];
  od;
  tgmat:=LeechBasis*ttgmat*LeechBasisinv;
  if not IsIntMat(tgmat) then beep(77112); fi;
  if TMTTmult(tgmat, GramLeech)<>GramLeech then beep(11241); fi;
  return(tgmat);
end;

#################3

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


####################33