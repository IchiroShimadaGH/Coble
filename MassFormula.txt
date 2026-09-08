#Read("MassFormula.txt");



ZpDiagonalize:=function(Gram, p)
  #
  local n,m,remainingG,newb,newbasis,newbasisup,newbasisdown,diags,
  ppower, nus,minnu,aa,ii,jj,nusmat,doneflag, eii, newsmallB, tG,
  kk, xy, ejj, AA, AAinv, zpdiagsrec, tqs, newdiags, tq, fq, xx, smallfqs,
  newdiagrec;
  #
  n:=Length(Gram);
  remainingG:=List(Gram);
  m:=n; #sizeremainingG;
  newbasis:=IdentityMat(n);
  diags:=[];
  ppower:=1;
  while m>0 do
    #Printn("before", remainingG, diags);
    #
    nus:=List(Flat(remainingG), aa->RatOrdpZ(aa, p));
    minnu:=Minimum(nus);
    if minnu>0 then
      ppower:=ppower*p^minnu;
      remainingG:=p^(-minnu)*remainingG;
      nus:=List(nus, nu->nu-minnu);
    fi;
    nusmat:=ChopVectToMat(nus, m, m);
    doneflag:=false;
    for ii in [1..m] do
      if nusmat[ii][ii]=0 then
        aa:=remainingG[ii][ii];
        eii:=MakeVectei(m, ii);
        newsmallB:=[eii];
        for jj in [1..m] do
          if jj<>ii then
            newb:=MakeVectei(m, jj)-remainingG[ii][jj]/aa*eii;
            Add(newsmallB, newb);
          fi;
        od;
        tG:=TMTTmult(newsmallB, remainingG);
        if m>1 then
          if not IsZeroMat(SubMatrix(tG, [1], [2..m])) then
            beep(42374113);
          fi;
          remainingG:=SubMatrix(tG, [2..m], [2..m]);
        else
          remainingG:=[];
        fi;
        Add(diags, [ppower, [[tG[1][1]]]]);
        newbasisup:=List([1..n-m], ii->newbasis[ii]);
        newbasisdown:=List([n-m+1..n], ii->newbasis[ii]);
        newbasis:=StackMats(newbasisup, newsmallB*newbasisdown);
        m:=m-1;
        doneflag:=true;

      fi;
      if doneflag then break; fi; #from for ii in [1..m] do
    od;
    if not doneflag then
      for ii in [1..m] do
        for jj in [ii+1..m] do
          if nusmat[ii][jj]=0 then
            AA:=SubMatrix(remainingG, [ii, jj], [ii, jj]);
            AAinv:=InverseMat(AA);
            eii:=MakeVectei(m, ii);
            ejj:=MakeVectei(m, jj);
            newsmallB:=[eii, ejj];
            for kk in [1..m] do
              if kk<>ii and kk<>jj then
                xy:=AAinv*[remainingG[ii][kk], remainingG[jj][kk]];
                newb:=MakeVectei(m, kk)-xy[1]*eii-xy[2]*ejj;
                Add(newsmallB, newb);
              fi;
            od;
            tG:=TMTTmult(newsmallB, remainingG);
            if m>2 then
              if not IsZeroMat(SubMatrix(tG, [1..2], [3..m])) then
                beep(423743);
              fi;
              remainingG:=SubMatrix(tG, [3..m], [3..m]);
            else
              remainingG:=[];
            fi;
            Add(diags, [ppower, SubMatrix(tG, [1..2], [1..2]) ]);
            newbasisup:=List([1..n-m], ii->newbasis[ii]);
            newbasisdown:=List([n-m+1..n], ii->newbasis[ii]);
            newbasis:=StackMats(newbasisup, newsmallB*newbasisdown);
            m:=m-2;
            doneflag:=true;
          fi;
          if doneflag then break; fi; #from for jj in [ii+1..m] do
        od;
        if doneflag then break; fi; #from for ii in [1..m] do
      od;
    fi;
    if not doneflag then beep(4265737133); fi;
    #Printn("after", remainingG, diags);
  od;
  #
  tqs:=List(diags, xx->xx[1]);
  newdiags:=[];
  for tq in Set(tqs) do
    smallfqs:=List(Positions(tqs, tq), ii->diags[ii][2]);
    fq:=DiagonalMats(smallfqs);
    newdiagrec:=rec(
      q:=tq,
      fq:=fq,
      smallfqs:=smallfqs
    );
    Add(newdiags, newdiagrec);
  od;
  #
  zpdiagsrec:=rec();
  zpdiagsrec.p:=p;
  zpdiagsrec.Gram:=Gram;
  zpdiagsrec.newbasis:=newbasis;
  zpdiagsrec.diags:=newdiags;
  return(zpdiagsrec);
end;


checkZpDiagonalize:=function(zpdiagsrec)
  local beep, newGram, xx, diag, M, nus, nu, p, smallfqs;
  #
  beep:=function(beepnumb)
    localbeep("checkZpDiagonalize", beepnumb); Error();
  end;
  newGram:=DiagonalMats(List(zpdiagsrec.diags, xx->xx.q*xx.fq));
  if TMTTmult(zpdiagsrec.newbasis, zpdiagsrec.Gram)<>newGram then
    beep(3612623);
  fi;
  p:=zpdiagsrec.p;
  if RatOrdp(Determinant(zpdiagsrec.newbasis), p)<>0 then
    beep(57534);
  fi;
  for diag in zpdiagsrec.diags do
    M:=diag.fq;
    smallfqs:=diag.smallfqs;
    if p=2 then
      Printn("smallfqs", smallfqs  mod 8);
    fi;
    if not IsSubset([1,2], Set(smallfqs, Length)) then beep(57334); fi;
    nus:=List(Flat(M), aa->RatOrdpZ(aa, p));
    if true in List(nus, nu->nu<0) then beep(4633); fi;
    if RatOrdp(Determinant(M), p)<>0 then
      beep(5757);
    fi;
  od;
  #
  return(true);
end;


##########################


TripleOdds:=function(tG)
  local tB, mm,newtB, newtG, aa, v1, new2tB, new2tG, v2,
  new3tG, new3tB;
  tB:=IdentityMat(3);
  mm:=1;
  newtB:=[tB[2]+mm*tB[1], tB[3]+mm*tB[1], tB[1]];
  newtG:=TMTTmult(newtB, tG);
  if newtG[1][1] mod 2 =0 then beep(54473); fi;
  aa:=newtG[1][1];
  v1:=newtB[1];
  new2tB:=[v1, newtB[2]-(newtG[1][2]/aa)*v1, newtB[3]-(newtG[1][3]/aa)*v1];
  new2tG:=TMTTmult(new2tB, tG);
  if new2tG[1][2]<>0 or new2tG[1][3]<>0 then beep(442333); fi;
  if new2tG[2][2] mod 2 =0 then beep(74274473746374637463); fi;
  aa:=new2tG[2][2];
  v1:=new2tB[1];
  v2:=new2tB[2];
  new3tB:=[v1, v2, new2tB[3]-(new2tG[2][3]/aa)*v2];
  new3tG:=TMTTmult(new3tB, tG);
  if not IsDiagonalMat(new3tG) then beep(447746); fi;
  if new3tG[3][3] mod 2 =0 then beep(47473); fi;
  new3tG:=new3tG mod 8;
  if (Determinant(new3tG)-Determinant(tG)) mod 8<>0  then beep(47473); fi;
  return(List([1,2,3], ii->[[new3tG[ii][ii]]]));
end;

_zpdrec2jdprec_two:=function(zpdrec)
  #
  local Jdqrec, p, nus, xx, minnu, maxnu,
  components, nu, jdqrec,  pos, diagrec,
  fq, detfq, d, smallfqs, aa, type, ones, twos, nn1, nn3, smallfq, esmallfq,
  tM, tripleodds, octane, lcomps, ii;
  #
  p:=2;
  Jdqrec:=rec();
  Jdqrec.p:=p;
  #
  nus:=List(zpdrec.diags, xx->RatOrdp(xx.q, p));
  if not IsDuplicateFreeList(nus) then beep(522323); fi;
  if CopySort(nus)<>nus then beep(5111272); fi;
  minnu:=nus[1];
  maxnu:=nus[Length(nus)];
  #
  components:=[];
  for nu in [minnu-1..maxnu+1] do
    jdqrec:=rec(nu:=nu, p:=p);
    Add(components, jdqrec);
  od;
  #
  for jdqrec in components do
    nu:=jdqrec.nu;
    if not nu in nus then
      #
      jdqrec.n:=0;
      jdqrec.type:="even";
      jdqrec.octane:=0;
      jdqrec.freebound:="not yet";
      #
    else
      pos:=SinglePosition(nus, nu);
      diagrec:=zpdrec.diags[pos];
      if diagrec.q<>p^nu then beep(61115433); fi;
      fq:=diagrec.fq;
      smallfqs:=diagrec.smallfqs;
      detfq:=Determinant(fq);
      if RatOrdp(detfq, p)<>0 then beep(15251818); fi;
      d:=detfq mod 8;
      #
      jdqrec.n:=Length(fq);
      #
      if 1 in List(smallfqs, Length) then
        type:="odd";
      else
        type:="even";
      fi;
      jdqrec.type:=type;
      #
      if jdqrec.type="even" then
        if d in [1,7] then octane:=0;
        elif d in [3,5] then octane:=4;
        else beep(57777);
        fi;
      elif jdqrec.type="odd" then
        ones:=[];
        twos:=[];
        for smallfq in smallfqs do
          if Length(smallfq)=1 then Add(ones, smallfq mod 8);
          elif Length(smallfq)=2 then  Add(twos, smallfq mod 8);
          else beep(1993658);
          fi;
        od;
        for esmallfq in twos do
          tM:=DiagonalMats([Remove(ones), esmallfq]);
          tripleodds:=TripleOdds(tM);
          Append(ones, tripleodds);
        od;
        nn1:=0; nn3:=0;
        for smallfq in ones do
          aa:=smallfq[1][1] mod 4;
          if aa=1 then nn1:=nn1+1;
          elif aa=3 then nn3:=nn3+1;
          else beep(746211);
          fi;
        od;
        octane:=(nn1-nn3) mod 8;
        #
      fi;
      #
      jdqrec.octane:=octane;
      jdqrec.freebound:="not yet";
    fi;
  od;
  #
  lcomps:=Length(components);
  if lcomps<3 then beep(4319135); fi;
  if components[2].type="odd" then
    components[1].freebound:="bound";
  else
    components[1].freebound:="free";
  fi;
  if components[lcomps-1].type="odd" then
    components[lcomps].freebound:="bound";
  else
    components[lcomps].freebound:="free";
  fi;
  for ii in [2..lcomps-1] do
    if components[ii-1].type="odd" or components[ii+1].type="odd" then
      components[ii].freebound:="bound";
    else
      components[ii].freebound:="free";
    fi;
  od;
  #
  Jdqrec.components:=components;
  return(Jdqrec);
end;

_zpdrec2jdprec_odd:=function(zpdrec)
  #
  local p, Jdqrec, nus, minnu, maxnu, components, nu,
  xx, jdqrec, pos, diagrec, fq, detfq;
  p:=zpdrec.p;
  Jdqrec:=rec();
  Jdqrec.p:=p;
  #
  nus:=List(zpdrec.diags, xx->RatOrdp(xx.q, p));
  if not IsDuplicateFreeList(nus) then beep(522323); fi;
  if CopySort(nus)<>nus then beep(5111272); fi;
  minnu:=nus[1];
  maxnu:=nus[Length(nus)];
  #
  components:=[];
  for nu in [minnu..maxnu] do
    jdqrec:=rec(nu:=nu, p:=p);
    Add(components, jdqrec);
  od;
  #
  for jdqrec in components do
    nu:=jdqrec.nu;
    if not nu in nus then
      #
      jdqrec.n:=0;
      jdqrec.d:=1;
      #
    else
      pos:=SinglePosition(nus, nu);
      diagrec:=zpdrec.diags[pos];
      if diagrec.q<>p^nu then beep(61115433); fi;
      fq:=diagrec.fq;
      jdqrec.n:=Length(fq);
      detfq:=Determinant(fq);
      if RatOrdp(detfq, p)<>0 then beep(15251818); fi;
      jdqrec.d:=detfq mod p;
    fi;
  od;
  #
  Jdqrec.components:=components;
  return(Jdqrec);
end;


ZpdrecTojdprec:=function(zpdrec)
  #
  local p;
  p:=zpdrec.p;
  if p=2 then return(_zpdrec2jdprec_two(zpdrec));
  else return(_zpdrec2jdprec_odd(zpdrec)); fi;
  beep(46463); return();
end;

######################


SpeaciesToMp:=function(p, species)
  #
  local ss, theMp, tt;
  #
  if species.sgn=0 then
    if species.val mod 2=0 then beep(553883);fi;
    ss:=(species.val+1)/2;
    theMp:=1/2;
    for tt in [1..ss-1] do
      theMp:=theMp/(1-p^(-2*tt));
    od;
  else
    ss:=species.val/2;
    if ss=0 then
      theMp:=1;
    else
      theMp:=1/2;
      for tt in [1..ss-1] do
        theMp:=theMp/(1-p^(-2*tt));
      od;
      theMp:=theMp/(1-species.sgn*p^(-ss));
    fi;
  fi;
  return(theMp);
end;

DiagonalMp:=function(jdqrec) #Mp
  local p, n, species, qres,
  freebound, octane, type, dim, t, d;
  #
  p:=jdqrec.p;
  n:=jdqrec.n;
  species:=rec();
  #
  if p=2 then
    freebound:=jdqrec.freebound;
    octane:=jdqrec.octane;
    type:=jdqrec.type;
    dim:=jdqrec.n;
    if type="even" then
      t:=dim/2;
    elif type="odd" then
      if dim mod 2 =0 then
        t:=(dim-2)/2;
      else
        t:=(dim-1)/2;
      fi;
    else
      beep(65334);
    fi;
    #
    if (octane in [0,1,7]) and freebound="free" then
      species.val:=2*t;
      species.sgn:=1;
    elif (octane in [2,6]) or freebound="bound" then
      species.val:=2*t+1;
      species.sgn:=0;
    elif (octane in [3,4,5]) and freebound="free" then
      if t=0 then beep(4366416656654165356656541654); fi;
      species.val:=2*t;
      species.sgn:=-1;
    else
      beep(5533);
    fi;
  else #p is odd
    d:=jdqrec.d;
    species:=rec();
    species.val:=n;
    if n mod 2=1 then
      species.sgn:=0;
    else
      qres:=Legendre(d mod p, p);
      if p mod 4=3 then
        if qres=1 then
          if n mod 4=0 then
            species.sgn:=1;
          elif n mod 4=2 then
            species.sgn:=-1;
          fi;
        elif qres=-1 then
          if n mod 4=0 then
            species.sgn:=-1;
          elif n mod 4=2 then
            species.sgn:=1;
          fi;
        fi;
      elif p mod 4=1 then
        if qres=1 then
          if n mod 4=0 then
            species.sgn:=1;
          elif n mod 4=2 then
            species.sgn:=1;
          fi;
        elif qres=-1 then
          if n mod 4=0 then
            species.sgn:=-1;
          elif n mod 4=2 then
            species.sgn:=-1;
          fi;
        fi;
      fi;
    fi;
  fi;
  #
  return(SpeaciesToMp(p, species));
end;


Massp:=function(Jdprec) #m_p(f)
  local p, diagonalproduct, jdqrec, diagonalfactor,
  thempf,  crossproductexp, ii, jj, nuii, nujj, nii, njj,
  jdqrecii, jdqrecjj, crossfactorexp, typefactor,
  nextjdqrec, neven, noddodd;
  #
  p:=Jdprec.p;
  #
  diagonalproduct:=1;
  for jdqrec in Jdprec.components do
    diagonalfactor:=DiagonalMp(jdqrec);
    diagonalproduct:=diagonalproduct*diagonalfactor;
  od;
  #
  crossproductexp:=0;
  for ii in [1..Length(Jdprec.components)] do
    jdqrecii:=Jdprec.components[ii];
    nuii:=jdqrecii.nu;
    nii:=jdqrecii.n;
    for jj in [ii+1..Length(Jdprec.components)] do
      jdqrecjj:=Jdprec.components[jj];
      nujj:=jdqrecjj.nu;
      njj:=jdqrecjj.n;
      #
      crossfactorexp:=(nujj-nuii)*(nii*njj)/2;
      crossproductexp:= crossproductexp+crossfactorexp;
      #
    od;
  od;
  #
  #
  thempf:=rec(ratfactor:=diagonalproduct, sqrootfactor:=p^(2*crossproductexp));
  if p=2 then
    neven:=0;
    noddodd:=0;
    for ii in [1..Length(Jdprec.components)] do
      jdqrecii:=Jdprec.components[ii];
      if jdqrecii.type="even" then neven:=neven+jdqrecii.n;
      else
        if ii<Length(Jdprec.components) then
          nextjdqrec:=Jdprec.components[ii+1];
          if nextjdqrec.type="odd" then
            noddodd:=noddodd+1;
          fi;
        fi;
      fi;
    od;
    typefactor:=2^(noddodd-neven);
    thempf.ratfactor:=thempf.ratfactor*typefactor;
  fi;
  #
  return(thempf);
end;

############################


TheJacobiD:=function(D)
  #
  local tJacobi, pp, vals, kkflag, dd, kk1flag, kk1, psi,
  flag, ii, kksptimes, kk, dd1, cdd, aa, ii2;
  #
  pp:=8*Product(SetMinus(Set(Factors(AbsInt(D))), [2]));;
  tJacobi:=function(mm)
    if mm mod 2 =0 then return(0);fi;
    return(Jacobi(D, mm));
  end;
  vals:=List([1..pp], tJacobi);
  #Printn(pp, vals);
  kkflag:=false;
  for dd in DivisorsInt(pp) do
    if dd<>pp then
      cdd:=pp/dd;
      flag:=true;
      for ii in [1..dd] do
        for aa in [1..cdd-1] do
          if vals[ii]<>vals[ii+aa*dd] then
            flag:=false; break; #from for ss in [1..cdd-1] do
          fi;
        od;
        if flag=false or (Gcd(ii, dd)>1 and vals[ii]<>0) then
          flag:=false; break; #from for ii in [1..dd] do
        fi;
      od;
      if flag then
        kk:=dd;
        kkflag:=true;
        break; #from for dd in DivisorsInt(pp) do
      fi;
    fi;
  od;
  if not kkflag then kk:=pp; fi;
  #
  kksptimes:=[];
  for ii in [1..kk] do
    if Gcd(ii, kk)=1 then
      Add(kksptimes, ii);
    fi;
  od;
  #
  kk1flag:=false;
  for dd1 in DivisorsInt(kk) do
    if dd1<>kk then
      flag:=true;
      for ii in kksptimes do
        #if (ii mod dd1)=1 and vals[ii]<>1 then
        if ((ii-1) mod dd1)=0 and vals[ii]<>1 then #for the case when dd1=1
          flag:=false;
          break; #from for ii in kksptimes do
        fi;
      od;
      if flag then
        kk1:=dd1;
        kk1flag:=true;
        break; #from for dd1 in DivisorsInt(kk1) do
      fi;
    fi;
  od;
  if not kk1flag then kk1:=kk; fi;
  #
  psi:=[];
  for ii in [1..kk1] do
    if Gcd(ii, kk1)<>1 then Add(psi, 0);
    else
      flag:=false;
      for ii2 in kksptimes do
        if (ii-ii2) mod kk1=0 then
          flag:=true;
          Add(psi, vals[ii2]);
          break; #from for ii2 in kksptimes do
        fi;
      od;
      if not flag then beep(53344); fi;
    fi;
  od;
  return(rec(D:=D, kk:=kk, kk1:=kk1, psi:=psi));
end;

checkTheJacobiD:=function(JacobiDrec, checkfactor)
  local D,M,jacobi,kk,kk1,psi,psilist,chi1, m;
  D:=JacobiDrec.D;
  M:=checkfactor*8*AbsInt(D);
  jacobi:=function(m)
    if m mod 2 =0 then return(0); fi;
    return(Jacobi(D, m));
  end;
  kk:=JacobiDrec.kk;
  kk1:=JacobiDrec.kk1;
  psilist:=JacobiDrec.psi;
  chi1:=function(m)
    if Gcd(m, kk)=1 then return(1); fi;
    return(0);
  end;
  psi:=function(m)
    return(psilist[((m-1) mod kk1)+1]);
  end;
  for m in [1..M] do
    if jacobi(m)<>chi1(m)*psi(m) then beep(4733); fi;
  od;
  Printn("checkTheJacobiD okay", D);
end;

####################


BernoulliVar:=X(Rationals, "x");

BernoulliPolynomial:=function(n)
  local k, f;
  f:=0*BernoulliVar;
  for k in [0..n] do
    f:=f+NrCombinations([1..n], k)*Bernoulli(k)*BernoulliVar^(n-k);
  od;
  return(f);
end;

TheZetaD:=function(D, s)
  #
  local JacobiD, kk, kk1, psilist, psi, Bs, zetaD,
  fac1, fac2,  fac3, p, I, aa, ee, r, psir, gausssq, gauss, gausssgn, fgauss,
  trigs, alpha, totalsgn, g;
  #
  JacobiD:=TheJacobiD(D);
  kk:=JacobiD.kk;
  kk1:=JacobiD.kk1;
  psilist:=JacobiD.psi;
  psi:=function(m)
    return(psilist[((m-1) mod kk1)+1]);
  end;
  Bs:=BernoulliPolynomial(s);
  #
  fac1:=1;
  for p in DivisorsInt(kk) do
    if IsPrime(p) then
      fac1:=fac1*(1-psi(p)/p^s);
    fi;
  od;
  #
  I:=E(4);
  aa:=I^(-s)+psi(-1)*I^s;
  if SignInt(aa^2)<>SignInt(D) then beep(4271945); fi;
  if aa=2 then alpha:=2;
  elif aa=-2 then alpha:=-2;
  elif aa=2*I then alpha:=2;
  elif aa=-2*I then alpha:=-2;
  else
    beep(447333121464);
  fi;
  fac2:=-2^s/(Factorial(s)*alpha);
  #
  fac3:=0;
  gauss:=0;
  fgauss:=.0;
  ee:=E(kk1);
  if SignInt(D)=1 then
    trigs:=List([1..kk1], jj->Cos(2*FLOAT.PI*jj/kk1));
  else
    trigs:=List([1..kk1], jj->Sin(2*FLOAT.PI*jj/kk1));
  fi;
  for r in [1..kk1] do
    psir:=psi(r);
    fac3:=fac3+psir*Value(Bs, [BernoulliVar], [r/kk1]);
    gauss:=gauss+psir*(ee^r);
    fgauss:=fgauss+psir*trigs[r];
  od;
  gausssq:=gauss^2;
  g:=AbsoluteValue(gausssq);
  if not IsInt(gausssq) then beep(43345); fi;
  if SignInt(gausssq)<>SignInt(D) then beep(42345); fi;
  if AbsoluteValue(AbsInt(gausssq)-fgauss^2)>.001 then beep(553344); fi;
  #
  if fgauss>.9 then gausssgn:=1;
  elif fgauss<-.9 then gausssgn:=-1;
  else beep(5473);
  fi;
  totalsgn:=gausssgn*SignInt(D);
  #
  zetaD:=rec(D:=D, s:=s, pifactor:=s, ratfactor:=fac1*fac2*fac3*totalsgn,
  g:=g, sqrootfactor:=1/g); #zetaD=pi^pifactor*ratfactor/sqrt(g);
  return(zetaD);
end;


###############################



DoubleFactorial:=function(m)
  local aa, xx;
  if m<0 then beep(77223); fi;
  aa:=1;
  for xx in [1..m] do
    if (xx-m) mod 2 =0 then
      aa:=aa*xx;
    fi;
  od;
  return(aa);
end;

TheGamma:=function(a)
  local beep, Grec, m, rat;
  #
  beep:=function(beepnumb)
    localbeep("theGamma", beepnumb); Error();
  end;
  #
  if a<1/2 then beep(2525); fi;
  if not IsInt(2*a) then beep(4763); fi;
  if IsInt(a) then
    Grec:=rec(ratfactor:=Factorial(a-1), pifactor:=0);
  elif IsInt(a-1/2) then
    m:=a-1/2;
    if m=0 then rat:=1;
    else
      rat:=DoubleFactorial(2*m-1)/2^m;
    fi;
    Grec:=rec(ratfactor:=rat, pifactor:=1/2);
  else
    beep(32632);
  fi;
  return(Grec);
end;

TheZeta:=function(m)
  local n, rat, zetarec;
  if m<=1 then beep(4573); fi;
  if m mod 2 <>0 then beep(5374); fi;
  n:=m/2;
  rat:=(-1)^(n+1)*Bernoulli(m)*2^m/2/Factorial(m);
  zetarec:=rec(pifactor:=m, ratfactor:=rat);
  return(zetarec);
end;


StdMass:=function(n, disc)
  #
  local s, D, gg, tt, zz, gammafac, jj;
  #
  if n mod 2 =0 then s:=n/2;
  else s:=(n+1)/2;
  fi;
  D:=(-1)^s*disc;
  gg:=rec(pifactor:=-n*(n+1)/4, ratfactor:=2, sqrootfactor:=1);
  #
  for jj in [1..n] do
    gammafac:=TheGamma(jj/2);
    gg.pifactor:=gg.pifactor+gammafac.pifactor;
    gg.ratfactor:=gg.ratfactor*gammafac.ratfactor;
  od;
  #
  for tt in [1..s-1] do
    zz:=TheZeta(2*tt);
    gg.pifactor:=gg.pifactor+zz.pifactor;
    gg.ratfactor:=gg.ratfactor*zz.ratfactor;
  od;
  #
  #Printn(n, gg.pifactor, gg.ratfactor); #Table 3
  #
  if n mod 2=0 then
    zz:=TheZetaD(D, s);
    gg.pifactor:=gg.pifactor+zz.pifactor;
    gg.ratfactor:=gg.ratfactor*zz.ratfactor;
    gg.sqrootfactor:=zz.sqrootfactor;
  fi;
  #
  if gg.pifactor<>0 then beep(424473); fi;
  #
  return(gg);
  #
end;


StdMassp:=function(n, d, p)
  local s, D, ee, stdp, tt;
  stdp:=1/2;
  if n mod 2 =0 then
    s:=n/2;
    D:=(-1)^s*d;
    if (2*d) mod p=0 then
      ee:=0;
    else
      ee:=Legendre(D, p);
    fi;
    stdp:=stdp/(1-ee*p^(-s));
  else
    s:=(n+1)/2;
  fi;
  for tt in [1..s-1] do
    stdp:=stdp/(1-p^(-2*tt));
  od;
  return(stdp);
end;



###################################



Mass:=function(tGram) #mf
  local n, d, ps, tstdmass, psfactor,p,
  zpdrec,Jdprec,tmpf,rat,sqrt,tstdmassp;
  n:=Length(tGram);
  d:=DeterminantIntMat(tGram);
  ps:=Set(Factors(2*d));
  ps:=SetMinus(ps, [1]);
  tstdmass:=StdMass(n, d);
  psfactor:=rec(ratfactor:=1, sqrootfactor:=1);
  for p in ps do
    zpdrec:=ZpDiagonalize(tGram, p);
    Jdprec:=ZpdrecTojdprec(zpdrec);
    tmpf:=Massp(Jdprec);
    tstdmassp:=StdMassp(n, d, p);
    psfactor.ratfactor:=psfactor.ratfactor*tmpf.ratfactor/tstdmassp;
    psfactor.sqrootfactor:=psfactor.sqrootfactor*tmpf.sqrootfactor;
  od;
  rat:=tstdmass.ratfactor*psfactor.ratfactor;
  sqrt:=SqRoot(tstdmass.sqrootfactor*psfactor.sqrootfactor);
  if not IsRat(sqrt) then beep(7463746367);fi;
  return(rat*sqrt);
end;

###################################
# Add 2022/02/14


RandomKneserpNbrL:=function(GramL, p)
  local dim, beep, cs, cc, findvv, findww, vv, ww,
  vvnrm, aa, bb, vvdual, snf, Lv, genvs, newGram;
  #
  beep:=function(beepnumb)
    localbeep("RandomKneserpNbrL", beepnumb); Error();
  end;
  #
  dim:=Length(GramL);
  cs:= [-1,0,1];
  cc:=1;
  findvv:=false;
  while not findvv  do
    vv:=RandomVectfromL(dim, cs);
    if IsZeroVect(vv mod p) then
      continue; #in while not findvv  do
    fi;
    vvnrm:=vv*GramL*vv;
    if vvnrm mod p=0 then
      aa:=vvnrm/p;
      if aa mod p=0 then
        findvv:=true;
        break; #from   while not findvv  do
      fi;
      findww:=false;
      while not findww  do
        ww:=RandomVectfromL(dim, cs);
        bb:=ww*GramL*vv;
        if bb mod p <>0 then
          findww:=true;
          cc:=-aa/2/bb mod p;
          vv:=vv+p*cc*ww;
          findvv:=true;
          break; #from while not findww  do
        fi;
      od;
    fi;
    cc:=cc+1;
    Append(cs, [-cc, cc]);
  od;
  vvnrm:=vv*GramL*vv;
  if IsZeroVect(vv mod p) then beep(6922); fi;
  if vvnrm mod p^2 <>0 then beep(698269823); fi;
  vvdual:=vv*GramL;
	snf:=SmithNormalFormIntegerMatTransforms(TransposedMat([vvdual]));
	if snf.normal[1][1] mod p =0  then beep(1112); fi;
  Lv:=snf.rowtrans;
  Lv[1]:=p*Lv[1];
  genvs:=p*Lv;
  Add(genvs, vv);
  genvs:=RemoveZeroVectsIntMat(HermiteNormalFormIntegerMat(genvs));
  if Rank(genvs)<>dim then beep(6986983); fi;
  newGram:=TMTTmult(genvs, GramL)/(p^2);
  if not IsEvenLattice(newGram) then beep(563773); fi;
  if DeterminantIntMat(newGram)<>DeterminantIntMat(GramL) then beep(9698263); fi;
  return(newGram);
end;
