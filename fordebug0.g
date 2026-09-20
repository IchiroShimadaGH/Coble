tGram:=[ [ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ],
 [ 2, -2, 2, 2, 0, 2, 2, 2, 2, 2, 0, 2 ], 
  [ 2, 2, -2, 0, 2, 2, 2, 2, 2, 4, 2, 4 ], [ 2, 2, 0, -2, 0, 2, 2, 2, 2, 1, 0, 1 ], 
  [ 2, 0, 2, 0, -2, 0, 2, 2, 2, 1, 2, 3 ], [ 2, 2, 2, 2, 0, -2, 0, 2, 2, 1, 4, 2 ], 
  [ 2, 2, 2, 2, 2, 0, -2, 0, 2, 0, 2, 0 ], [ 2, 2, 2, 2, 2, 2, 0, -2, 0, 2, 2, 2 ], 
  [ 2, 2, 2, 2, 2, 2, 2, 0, -2, 4, 0, 0 ], [ 2, 2, 4, 1, 1, 1, 0, 2, 4, -2, 4, 4 ], 
  [ 2, 0, 2, 0, 2, 4, 2, 2, 0, 2, -2, 1 ], [ 2, 2, 4, 1, 3, 2, 0, 2, 0, 1, 2, -2 ] ];



  
  tModLatticeKerRec:=function(bigGram)
	#
	local  beep, bigleng, id, hn, rk, basis, basisKer, ii, redGram, T,
  trec, phi, zzredGram, diff;
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
  diff:=bigleng-rk;
	hn:=HermiteNormalFormIntegerMatTransform(bigGram);
	T:=hn.rowtrans;
	basis:=List([1..rk], ii->T[ii]);
	basisKer:=List([rk+1..bigleng], ii->T[ii]);
	redGram:=TMTTmult(basis, bigGram);
  zzredGram:=MatMatToMat(
      [[redGram, NullMat(rk, diff)], 
        [NullMat(diff, rk), NullMat(diff, diff)]]
      );
  if TMTTmult(T, bigGram)<>zzredGram then beep(1481472); fi;
  if CokerTorsion(basis)<>[] then beep(613321); fi;
  if basisKer=[] then beep(333321); fi;
  if CokerTorsion(basisKer)<>[] then beep(615221); fi;
  if not IsZeroMat(basisKer*bigGram) then beep(225221); fi;
  phi:=SubMatrix(InverseMat(T), [1..bigleng], [1..rk]);
  if TMTTmult(phi, redGram)<>bigGram then beep(616996); fi;
  trec:=rec(
    redGram:=redGram, 
    basis:= basis,
    basisKer:=basisKer, 
    phi:=phi
  );
	return(trec);
end;


  tModLatticeKerRec(tGram);


  ######