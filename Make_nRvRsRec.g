#Read("Make_nRvRsRec.g");



Make_nRvRsRec:=function(GramR, embR, GramL)
	#
	# make the record for vectors in GramRdual
  # with -2<=squarenorm <=0
  # and their embeddings in L26
	#
	local beep, svrec, nRs, vRss, vv, vvnR,  nRvRsrec, pos, emb2dual,
  vRInLss, vRs,  GramRdual, LGrec, dd;
  #
  beep:=function(beepnumb)
    localbeep("Make_nRvRsRec", beepnumb); Error();
  end;
  #
  GramRdual:=InverseMat(GramR); #GramR is defined outside 
  emb2dual:=GramRdual*embR;
  dd:=Llcm(List(Flat(GramRdual), DenominatorRat));
  #
	svrec:=ShortestVectors((-dd)*GramRdual, 2*dd);
	if svrec.vectors=[] then beep(72727272); fi;
  #
  nRs:=(-1/dd)*Set(svrec.norms);
  Add(nRs, 0);
  Sort(nRs);
	#
  vRss:=List(nRs, nR->[]);
  Add(vRss[SinglePosition(nRs, 0)], MakeZeroVect(Length(GramR)));
  for vv in svrec.vectors do
    vvnR:=vv*GramRdual*vv;
    #if vvnR=-2 then continue; fi;
    # We need the elenets of sqnorm -2 in calculating the adjcham.
    # Hence we delete if vvnR=-2 then continue; fi;
    pos:=SinglePosition(nRs, vvnR);
    Append(vRss[pos], [vv, (-1)*vv]);
  od;
  #
  vRInLss:=List(vRss, vRs->vRs*emb2dual);
  #
	#for check
	#
  if List(vRss, vRs->Set(vRs, vv->vv*GramRdual*vv))<>List(nRs, xx->[xx]) then
    beep(52527);
  fi;
  if false in List(vRss, IsDuplicateFree) then beep(16191); fi;
  if false in List(vRInLss, IsDuplicateFree) then beep(126191); fi;
  if List(vRInLss, vRs->Set(vRs, vv->vv*GramL*vv))<>List(nRs, xx->[xx]) then
    beep(52527);
  fi;
  #
  nRvRsrec:=rec(
    nRs:=nRs,
    vRss:=vRss,
    vRInLss:=vRInLss
  );
	return(nRvRsrec);
end; #the end of Make_nRvRsRec: