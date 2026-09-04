#Read("Make_nSaSvSss.g");

Make_nSaSvSss:=function(Borrec, nowweyl)
	#
	# nSaSvSss is a list of records
	# nSaSvSsrec:=rec(nS, aS, vSliftss).
	# vSliftss is a list of [vS, lifts].
	#
  local beep, GramL, GramS, GramSdual, dd, intGramSdual, nowweyldual,
  wS, wSdual, embSdual, nSaSvSss, Rpos, nR, nS,
  aRvRInLss, aR, vRInLs, aS, lifts, task, lift, vSliftss, vSlifts, vS,
  nSaSvSsrec, aRvRInLs, nSaSs, xx, LGrec, nRvRsrec, totalvRInLs;
  #
  beep:=function(beepnumb)
    localbeep("Make_nSaSvSss", beepnumb); Error();
  end;
  #
  nRvRsrec:=Borrec.nRvRsrec; # will be made in BorchardsL26AF
  GramL:=Borrec.GramL;
  GramS:=Borrec.GramS;
  GramSdual:=Borrec.GramSdual;
  dd:=Llcm(List(Flat(GramSdual), DenominatorRat));
  intGramSdual:=dd*GramSdual;
  if not IsIntMat(intGramSdual) then beep(88811); fi;
  nowweyldual:=nowweyl*GramL;
  wS:=nowweyl*Borrec.projS;
  if wS*GramS*wS<=0 then beep(619611); fi;
  wSdual:=nowweyl*Borrec.projSdual;
  if wSdual<>wS*GramS then beep(158756781); fi;
  if not IsIntVect(wSdual) then beep(58111); fi;
  #
  embSdual:=Borrec.embSdual;
	#in terms of the STANDARD basis of L and the dual basis of Sdual
  #
  nSaSvSss:=[];
  #
  for Rpos in [1..Length(nRvRsrec.nRs)] do
    nR:=nRvRsrec.nRs[Rpos];
    if nR=-2 then continue; fi;
    nS:=-2-nR;
    if not nS<0 then beep(69816981); fi;
    totalvRInLs:=nRvRsrec.vRInLss[Rpos];
    aRvRInLss:=NewClassifyByTags(totalvRInLs, vv->vv*nowweyldual);#in gaptools
    #
    for aRvRInLs in aRvRInLss do
      aR:=aRvRInLs[1];
      vRInLs:=aRvRInLs[2];
      aS:=1-aR;
      #
      lifts:=[];
      #
      task:=function(vInSdual)
        local vSInL, vRInL, vL;
        vSInL:=vInSdual*embSdual;
        for vRInL in vRInLs do
          vL:=vSInL+vRInL;
          if IsIntVect(vL) then
            Add(lifts, vL);
          fi;
        od;
      end;
      #
      AffES(intGramSdual, wSdual, dd*aS, dd*nS,  true, task);
      vSliftss:=NewClassifyByTags(lifts, vv->vv*Borrec.projS);##in gaptools
      for vSlifts in vSliftss do # for check
        vS:=vSlifts[1];
        if vS*GramS*vS<>nS then beep(62191); fi;
        if vS*wSdual<>aS then beep(62291); fi;
        for lift in vSlifts[2] do
          if lift*nowweyldual<>1 then beep(261112); fi;
          if lift*GramL*lift<>-2 then beep(91277); fi;
        od;
      od;
      #
      if vSliftss<>[] then
        nSaSvSsrec:=rec(
          nS:=nS,
          aS:=aS,
          vSliftss:=vSliftss
        );
        Add(nSaSvSss, nSaSvSsrec);
        #Printn("nS", nS, "aS", aS, "nops", Length(vSliftss));
      else 
        #Printn("nS", nS, "aS", aS, "nops", 0);
      fi;
    od;
  od; #the end of for Rpos in [1..Length(nRvRsrec.nRs)] do
  #
  nSaSs:=List(nSaSvSss, xx->[xx.nS, xx.aS]);
  SortParallel(nSaSs, nSaSvSss);
  return(nSaSvSss);
end; #the end of Make_nSaSvSss.


check_nSaSvSss:=function(Borrec,  tweyl, nSaSvSss)
	local beep, twdual, twS, navS, nS, aS, navSrec,
	vSlifts, vS, lifts,  vSdual, rL, vSliftss, twSdual, GramL;
	#
	beep:=function(beepn) 
    localbeep("check_nSaSvSss", beepn);Error(); 
  end;
	#
  GramL:=Borrec.GramL;
	twdual:=tweyl*GramL;
	twS:=tweyl*Borrec.projS;
  twSdual:=twS*Borrec.GramS;
	for navSrec in nSaSvSss do
		nS:=navSrec.nS;
		aS:=navSrec.aS;
    vSliftss:=navSrec.vSliftss;
		if nS>=0 then beep(2); fi;
		for vSlifts in vSliftss do
      vS:=vSlifts[1];
      lifts:=vSlifts[2];
      vSdual:=vS*Borrec.GramS;
      if vS*twSdual<>aS then beep(91919); fi;
      if vSdual*vS<>nS then beep(99991); fi;
      for rL in lifts do 
        if rL*Borrec.projS<>vS then beep(2311); fi;
        if rL*GramL*rL<>-2 then beep(5); fi;
        if rL*twdual<>1 then beep(221); fi; 
      od;
		od;
	od;
	#  
	return(true);
end;


ChamNumData:=function(nSaSvSss)
  local trec, nums;
  nums:=List(nSaSvSss, trec->[trec.nS, trec.aS, Length(trec.vSliftss)]);
  return(Collected(nums));
end;

####