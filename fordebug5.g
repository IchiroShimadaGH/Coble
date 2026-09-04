

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
        if rL*twdual<>1 then beep(221); fi; koko(1);
      od;
		od;
	od;
	#  
	return(true);
end;

check_nSaSvSss(Borrec, weyl0, nSaSvSss0);