#Read("Task20260829Ub.g");

#########################
Read("paperaffineConway/affineConwayCompdata.txt");
Read("OGLat.g");
Read("ImKerOLtoOqL.g");
Read("checktools.g");

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
    if cc=1000 then beep(182871); fi;
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


###
