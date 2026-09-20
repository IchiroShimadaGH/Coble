
#Read("SplitConsTools.g");


NpriFindRatssOnK3:=function(GramSX, ample, uptodeg)
  local rats, extradeg, ratsdual, task, tdeg, ratss, drats;
  #
  rats:=[];
  ratsdual:=[];
  tdeg:=0;
  ratss:=[];
  drats:=[];
  #
  task:=function(tv)
    local tvdual;
    for tvdual in ratsdual do
      if tvdual*tv<0 then return(); fi;
    od;
    Add(drats, tv);
    Add(ratsdual, tv*GramSX);
    return();
  end;
  #
  while tdeg<uptodeg do
    tdeg:=tdeg+1;
    drats:=[];
    AffES(GramSX, ample, tdeg, -2, true, task);
    Add(ratss, drats);
  od;
  #
  return(ratss);
end;



Cover22:=function(M)
  local mm, tv, ttv, tt, ii;
  mm:=Length(M);
  tv:=[List([1..mm], ii->2)];
  ttv:=TransposedMat(tv);
  tt:=[[2]];
  return(MatMatToMat( [[tt, tv], [ttv, M]]));
end;


ModLatticeKerRec:=function(bigGram)
	#
	local  beep, bigleng, id, hn, rk, basis, basisKer, ii, redGram, T,
  trec, phi;
	#
  beep:=function(beepnumb)
    localbeep("ModLatticeKerRec", beepnumb); Error();
  end;
	#
	bigleng:=Length(bigGram);
  rk:=Rank(bigGram);
  if rk=bigleng then 
    id:=IdentityMat(bigleng);
    trec:=rec(
      redGram:=CopyMat(bigGram), 
      basis:= id, 
      basisKer:=[], 
      phi:=id
    );
	  return(trec);
  fi;
  #
	hn:=HermiteNormalFormIntegerMatTransform(bigGram);
	T:=hn.rowtrans;
	basis:=List([1..rk], ii->T[ii]);
	basisKer:=List([rk+1..bigleng], ii->T[ii]);
	redGram:=TMTTmult(basis, bigGram);
  if CokerTorsion(basis)<>[] then beep(613321); fi;
  if basisKer=[] then beep(333321); fi;
  if CokerTorsion(basisKer)<>[] then beep(615221); fi;
  if not IsZeroMat(basisKer*bigGram) then beep(225221); fi;
  phi:=SubMatrix(InverseMat(T), [1..bigleng], [1..rk]);
  trec:=rec(
    redGram:=redGram, 
    basis:= basis,
    basisKer:=basisKer, 
    phi:=phi
  );
	return(trec);
end;


WGraphToGramh:=function(w, kk)
  #
  local ii, jj, tv, tGram, redGrec, redGram,
  th, thdual, gensplcons, M, trec, xx;
  #
  if IsList(w[1]) then 
    M:=List(w, tv->List(tv));
    for ii in [1..kk] do 
      M[ii][ii]:=-2;
    od;
  else  
    M:=(-2)*IdentityMat(kk);
    ii:=1; jj:=2;
    for xx in w do 
      M[ii][jj]:=xx;
      M[jj][ii]:=xx;
      jj:=jj+1;
      if jj>kk then 
        ii:=ii+1;
        jj:=ii+1;
      fi;
    od;
  fi;
  #
  tGram:=Cover22(M);
  redGrec:=ModLatticeKerRec(tGram);
  redGram:=redGrec.redGram;
  th:=MakeVectei(kk+1, 1)*redGrec.phi;
  thdual:=th*redGram;
  if th*thdual<>2 then beep(18821); fi;
  #
  #
  gensplcons:=[];
  for jj in [1..kk] do 
    tv:=MakeVectei(kk+1, jj+1)*redGrec.phi;
    if tv*redGram*tv<>-2 then beep(88111); fi;
    if tv*thdual<>2 then beep(281834); fi;
    Add(gensplcons, tv);
  od;
  if TMTTmult(gensplcons, redGram)<>M then beep(916621); fi;
  #
  #
  trec:=rec(
    kk:=kk, 
    adj:=M, 
    h:=th, 
    gensplcons:=gensplcons, 
    Gram:=redGram
  );
  return(trec);
end;

NoOverlattice:=function(Gram, h)
  local discrec, ii, discgvs, Gramdual,  discf, tv, ttv, newvs, H, newGram, newh;
  #
  Gramdual:=InverseMat(Gram);
  discrec:=DiscriminantForm(Gram);
  discgvs:=Cartesian(List(discrec.discg, ii->[0..ii-1]));
  discf:=discrec.discf;
  for tv in discgvs do 
    if IsZeroVect(tv) then continue; fi;
    if modtZ(tv*discf*tv)=0 then 
      ttv:=tv*discrec.reps_dual;
      newvs:=CopyAdd(Gram, ttv);
      H:=HermiteNormalFormIntegerMat(newvs);
      H:=Filtered(H, xx->not IsZeroVect(xx));
      newGram:=TMTTmult(H, Gramdual);
      newh:=h*Gram*InverseMat(H);
      if AffESstd(newGram, newh, 0, -2, true)=[] then 
        return([false, ttv*Gramdual]);
        #
        # We can add ttv*Gramdual 
        # without destroying the ampleness of h.
        #
      fi;
    fi;
  od;
  return(true);
end;



IsGeom:=function(tGram, th)
  local nn;
  nn:=Length(tGram);
  if Rank(tGram)<>nn then return([false, 1]); fi;
  if SignatureQ(tGram)<>[nn, 1, nn-1] then return([false, 2]); fi;
  if th*tGram*th<>2 then return([false, 3]); fi;
  if not PrimitivelyEmbeddableInK3Lattice(tGram) then return([false, 4]); fi;
  if AffESstd(tGram, th, 0, -2, true)<>[] then return([false, 5]); fi;
  if AffESstd(tGram, th, 1, 0, true)<>[] then return([false, 6]); fi;
  return(true);
end;


GetTotalSplcons:=function(Gram, h)
  #
  local tinvol, ratss, scons, dones, tr, ttr;
  #
  tinvol:=RatInvolRecP2(Gram, h, h).invol;
  ratss:= NpriFindRatssOnK3(Gram, h, 2);
  if ratss[1]<>[] then Error("there is a tritangent line"); buzz(47612); fi;
  scons:=[];
  dones:=[];
  for tr in ratss[2] do
    if tr in dones then continue; fi;
    ttr:=tr*tinvol;
    if ttr<>tr then 
      Add(scons, [tr, ttr]);
      Append(dones, [tr, ttr]);
    fi; 
  od;
  return(scons);
end;


  ####################