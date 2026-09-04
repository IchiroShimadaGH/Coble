#Read("Adj.g");


IsWeylVect:=function(Borrec, tweyl) # used in check_AdjChamsdata
	local mG, GramL;
  GramL:=Borrec.GramL;
	if not IsPrimitiveVect(tweyl) then return(false); fi;
	if not tweyl*GramL*tweyl=0 then return(false); fi;
	mG:=(-1)*(HoloSphereLattice(GramL, tweyl).Gram);
	if Length(ShortestVectors(mG, 2).vectors)>0 then return(false); fi;
	return(true);
end;

IsWalldefv:=function(Borrec, defv) 
  defvdual:=defv*Borrec.GramS;
  if defvdual*defv>=0 then beep(8888111); fi;
  if not IsIntVect(defvdual) then beep(17721); fi;
  if not IsPrimitiveVect(defvdual) then beep(22721); fi;
  return(true);
end;


IsMinusTwoWall:=function(Borrec, defv) 
  local sqnormtv, mmdefv, mm;
  sqnormtv:=defv*(Borrec.GramS)*defv;
	mm:=SqRoot(-2/sqnormtv);
  if not IsInt(mm) then return(false); fi;
  mmdefv:=mm*defv;
  if not IsIntVect(mmdefv) then return(false); fi;
  return(mmdefv);
end;

AdjWeyl:=function(Borrec, tweyl, defv)
  #
  local beep, GramL, nrms, mm, rLs, aa, aadefvL, aanrmv, nRvRsrec,
  pos, nR, ttv, tr,  rLsdual, rLsdual_divided_by_tw, newweyl, adjrec;
  #
  beep:=function(beepnumb)
		localbeep("AdjWeyl", beepnumb); Error();
	end;
	#
  GramL:=Borrec.GramL;
  defvdual:=defv*(Borrec.GramS);
  nrmv:=defvdual*defv;
  mm:=SqRoot(-2/nrmv);
  #
  # We make a list, modulo pm 1, of the roots r in L 
  # whose hyperplane contains the wall in question (the wall defined by defv in S). 
  # Since defvdual is primitive in Sdual, the dual vector r^\vee
  # can be written in the form r^\vee = aa*defv + v,
  # where aa is a nonnegative integer and v is an element of R^\vee
  # with norm between -2 and 0, inclusive.
  #
  rLs:=[];
  for aa in [0..mm] do
    aadefvL:=aa*defv*Borrec.embS;
    aanrmv:=aa*aa*nrmv;
    #
    #RecNames(Borrec.nRvRsrec); [ "nRs", "vRInLss", "vRss" ]
    #
    nRvRsrec:=Borrec.nRvRsrec;
    pos:=0;
    for nR in nRvRsrec.nRs do 
      pos:=pos+1;
      if nR+aanrmv=-2 then
        for ttv in nRvRsrec.vRInLss[pos] do
          tr:=aadefvL+ttv;
          if IsIntVect(tr) then
            if tr*GramL*tr<>-2 then beep(524132); fi;
            if IniValnz(tr)>0 then Add(rLs, tr); else Add(rLs, -tr); fi;
          fi;
        od;
      fi;
    od;
  od;
  #
  # We order the hyperplanes defined by the minus-two vectors 
  # in rLs according to the order in which they intersect a general line passing through tweyl.
  # For a sufficiently small positive real number a, 
  # we consider the line passing through tweyl and va, given by
  # -t tweyl + va,
  # where 
  # va:= [a, a^2, a^3, \dots, a^26].
  # The intersection with the hyperplane defined by <x, r>=0 is given by
  # t=  <va, r>/<w, r>
  # Therefore, from the particular form of va, 
  # it suffices to order the vectors obtained by dividing the dual-vector representation of r
  # by the inner product of w and r, using the lexicographic order.
  #
  rLs:=Set(rLs);
  rLsdual:=List(rLs, xr->xr*GramL);
  rLsdual_divided_by_tw:=List(rLsdual, xr->xr/(tweyl*xr));
  #
  if not IsDuplicateFreeList(rLsdual_divided_by_tw) then
    beep(661231);
  fi;
  #
  newweyl:=List(tweyl);
  for tr in rLs do
    newweyl:=newweyl+(newweyl*GramL*tr)*tr;
  od;
  #
  adjrec:=rec(
    defv:=List(defv),
    rLs:=rLs, 
    oldweyl:=List(tweyl);
    newweyl:=newweyl
  );
  return(adjrec);
end;


CheckAdjrec:=function(Borrec, adjrec, trialnumb)
  local beep, GramL, rLs, tr, tvinL, inttvinL, ttwall,
  orthowall, orthoroots, rankL, trial, ttrLs, ttrLsdual, fflag,
  u, randscs, xr, parats, ttnewweyl;
  #
  beep:=function(beepnumb)
		localbeep("CheckAdjrec", beepnumb); Error();
	end;
	#
  GramL:=Borrec.GramL;
  IsWalldefv(Borrec, adjrec.defv);
  IsWeylVect(Borrec, adjrec.newweyl);
  #
  rLs:=adjrec.rLs;
  for tr in rLs do
    if not IsIntVect(tr)  then beep(22511); fi;
    if tr*GramL*tr<>-2 then beep(42511); fi;
    if Rank([adjrec.defv, tr*Borrec.projS])<>1 then beep(33211); fi;
  od;
  #
  tvinL:=adjrec.defv*Borrec.embS;
  inttvinL:=Llcm(List(tvinL, DenominatorRat))*tvinL;
  ttwall:=OrthogonalCompRec(GramL, Concatenation([inttvinL], Borrec.embR));
  orthowall:=OrthogonalCompRec(GramL, ttwall.basis);
  if Length(orthowall.basis)<>Length(Borrec.GramR)+1 then beep(262367); fi;
  orthoroots:=ShortestVectors((-1)*orthowall.Gram, 2).vectors;
  orthoroots:=orthoroots*orthowall.basis;
  if Set(orthoroots, ttv->ttv*GramL*ttv)<>[-2] then beep(71726); fi;
  orthoroots:=Union(orthoroots, (-1)*orthoroots);
  if orthoroots<>Union(rLs, (-1)*rLs)) then beep(42211); fi;
  #
  rankL:=26;
  for trial in [1..trialnumb] do
    ttrLs:=List(rLs);
    ttrLs:=Shuffle(ttrLs);
    ttrLsdual:=ttrLs*GramL;
    fflag:=false;
    u:=RandomVect(rankL);
    while not fflag do
      u:=2*u+RandomVect(rankL);
      randscs:=List(ttrLsdual, xr->(u*xr)/(adjrec.oldweyl*xr));
      if IsDuplicateFreeList(randscs) then
        fflag:=true;
      fi;
    od;
    parats:=List(ttrLs, xr->(u*GramL*xr)/(tweyl*GramL*xr));
    SortParallel(parats, ttrLs);
    ttnewweyl:=List(adjrec.oldweyl);
    for tr in ttrLs do
      ttnewweyl:=ttnewweyl+(ttnewweyl*GramL*tr)*tr;
    od;
    if ttnewweyl<>adjrec.newweyl then
      beep(33221);
    fi;
  od;
  #
  return(true);
  #
end;

check_AdjChamsdata:=function(Borrec, oldweyl, torbits, wallsresult, adjdata)
	#
	local beep, ii, tr, newweyl, GramS, GramL;
	#
	beep:=function(beepn) 
    localbeep("wrong in check_AdjChamsdata", beepn); Error(); 
  end;
  #
	if Length(adjdata)<>Length(torbits) then
		beep(0);
	fi;
  #
  GramL:=Borrec.GramL;
  GramS:=Borrec.GramS;
  #
	for ii in [1..Length(adjdata)] do
		if adjdata[ii][1]=_nw then
			if wallsresult[ii][1]<>false then beep(1); fi;
		elif adjdata[ii][1]=_mtc then
			tr:=adjdata[ii][2];
			if tr*GramS*tr<>-2 then beep(2); fi;
			if tr*GramS*ampleX<=0 then beep(3); fi;
		elif adjdata[ii][1]=_newcham then
			newweyl:=adjdata[ii][2];
      if not IsWeylVect(GramL, newweyl) then beep(4);fi;
		fi;
	od;
	#
	return(true);
end;


Make_AdjChamsdata:=function(tweyl, torbits, wallsresult)
	local beep, adjdata, ii, tvdual, tv, mtcflag, sqnormtv, n,
	Rtvs, a, atvdualL, tnorm, nRvRinL,
	ttv, tr, Rtvsdual, fflag, u, newweyl, reflects, enteringwall,
	parats, ttRtvs, trial, ttnewweyl, randscs, Rtvsdual_divided_by_tw, xr,
	#midweyl, midweylproj, midflag,GramSdual,
  tvinL, inttvinL, orthowall, orthoroots, ttwall, GramS, GramL;
	#
	# midpt of weyl and weyl*opptg is not necessarily on SX!
	# Consider the sphereical A2! 2017/10/28
	#
	beep:=function(beepnumb)
		localbeep("wrong Make_AdjChamsdata", beepnumb); Error();
	end;
	#
  GramL:=Borrec.GramL;
  GramS:=Borrec.GramS;
  GramSdual:=Borrec.GramSdual;
  #
	adjdata:=[];
	for ii in [1..Length(wallsresult)] do
		tvdual:=torbits[ii][1];
    if not IsPrimitiveVect(tvdual) then beep(61521); fi;
		if not wallsresult[ii][1] then #not a wall
			Add(adjdata, [_nw, tvdual]);
		else
			tv:=tvdual*GramSdual;
			mtcflag:=false;
			sqnormtv:=tv*GramS*tv;
			n:=RootInt(FloorRat((-2)/sqnormtv));
			if n*n*sqnormtv=-2 then
				if IsIntVect(n*tv) then
					mtcflag:=true;
					Add(adjdata, [_mtc, n*tv]);
				fi;
			fi;
      #
      if mtcflag then continue; fi;
      # in for ii in [1..Length(wallsresult)] do
      #
      Rtvs:=[];
      #
      # We make a list, modulo \pm 1, of the roots r in L 
      # whose hyperplane contains the wall in question (the wall defined by tvdual in S). Since tvdual is primitive, the dual vector r^\vee
      # can be written in the form r^\vee = a*tvdual + v,
      # where a is a nonnegative integer and v is an element of R^\vee
      #  with norm between -2 and 0, inclusive.
      #
      for a in [0..n] do
        atvdualL:=a*tv*Borrec.embS;
        tnorm:=a*a*sqnormtv;
        for nRvRinL in thenRvRinLs do
          if nRvRinL[1]+tnorm=-2 then
            for ttv in nRvRinL[2] do
              tr:=atvdualL+ttv;
              if IsIntVect(tr) then
                #if IniValnz(tr)>0 then
                #if not -tr in Rtvs then # 2017 Oct 26 correct
                #
                # The norm of tr should be -2:
                #
                if IniValnz(tr)>0 then
                  Add(Rtvs, tr);
                else
                  Add(Rtvs, -tr);
                fi;
              fi;
            od;
          fi;
        od;
      od;
      Rtvs:=Set(Rtvs);
      Rtvsdual:=List(Rtvs, xr->xr*GramL);
      Rtvsdual_divided_by_tw:=List(Rtvsdual, xr->xr/(tweyl*xr));
      SortParallel(Rtvsdual_divided_by_tw, Rtvs);
      #
      if not IsDuplicateFreeList(Rtvsdual_divided_by_tw) then
        Printn("not DuplicateFreeList in Make_AdjChamsdata 1");
        Error();
      fi;
      #
      #
      #################
      #
      #Added 26 Oct 2017
      #
      #
      if checklevel>0 then
        #Printn("Added 2017/10/28");
        tvinL:=tv*Lrec.emb1;
        inttvinL:=Lcm(List(tvinL, DenominatorRat))*tvinL;
        ttwall:=OrthogonalComp(Concatenation([inttvinL], Lrec.emb2), GramL);
        orthowall:=OrthogonalComp(ttwall.basis, GramL);
        if Length(orthowall.basis)<>rankR+1 then beep(262367); fi;
        orthoroots:=ShortestVectors((-1)*orthowall.Gram, 2).vectors;
        orthoroots:=orthoroots*orthowall.basis;
        if Set(orthoroots, ttv->ttv*GramL*ttv)<>[-2] then beep(71726); fi;
        orthoroots:=Concatenation(orthoroots, (-1)*orthoroots);
        if not IsEqualSet(orthoroots, Concatenation(Rtvs, (-1)*Rtvs)) then
          Printn("Rtvs is wrong in Make_AdjChamsdata");
          Error();
        fi;
      fi;
      #
      ################
      #
      #
      #fflag:=false;
      #while not fflag do
      #	u:=RandomVect(rankL);
      #	randscs:=List(Rtvsdual, xr->(u*xr)/(tweyl*xr));
      #	if IsDuplicateFreeList(randscs) then
      #		fflag:=true;
      #	fi;
      #od;
      #
      # the random vect u above will not be used: 2017/10/05
      #
      ttRtvs:=CopyMat(Rtvs); #this will be used in the check below:
      #
      #parats:=List(Rtvs, xr->(u*GramL*xr)/(tweyl*GramL*xr));
      #if not IsDuplicateFreeList(parats) then
      #	Printn("not DuplicateFreeList in Make_AdjChamsdata");
      #	Error();
      #fi;
      #
      #SortParallel(parats, Rtvs);
      # parats will not be used: 2017/10/05
      #
      SortParallel(Rtvsdual_divided_by_tw, Rtvs);
      #Be careful about the order of the args.
      #We want to sort Rtvs ACCORDING to parats.
      #
      newweyl:=CopyVect(tweyl);
      #
      #Printn("Rtvs", Length(Rtvs), ii);
      for tr in Rtvs do
        if Rank([tr*Lrec.proj1, tv])<>1 then
          Printn("tr should be perp to the wall in Make_AdjChamsdata");
          Error();
        fi;
        if tr*GramL*tr<>-2 then
          Printn("tr should be -2 in Make_AdjChamsdata");
          Error();
        fi;
        newweyl:=newweyl+(newweyl*GramL*tr)*tr;
      od;
      #
      enteringwall:=-tvdual;
      #
      Add(adjdata, [_newcham, newweyl, enteringwall]);
      #
      if checklevel>0 then
        for trial in [1..adjchecknumb] do
          fflag:=false;
          u:=RandomVect(rankL);
          while not fflag do
            u:=2*u+RandomVect(rankL);
            randscs:=List(Rtvsdual, xr->(u*xr)/(tweyl*xr));
            if IsDuplicateFreeList(randscs) then
              fflag:=true;
            fi;
          od;
          parats:=List(ttRtvs, xr->(u*GramL*xr)/(tweyl*GramL*xr));
          SortParallel(parats, ttRtvs);
          ttnewweyl:=CopyVect(tweyl);
          for tr in ttRtvs do
            ttnewweyl:=ttnewweyl+(ttnewweyl*GramL*tr)*tr;
          od;
          if ttnewweyl<>newweyl then
            Printn("Different newweyls in Make_AdjChamsdata");
            Error();
          else
            #Printn("newweyl ok: added  2017 Oct 26"); # added  2017 Oct 26
          fi;
        od;
        if adjchecknumb>0 then
          #Printn("newweyl ok done: ii=", ii, "adjchecknumb=", adjchecknumb); # added  2017 Oct 26
        fi;
      fi; #if checklevel>0 then
      #
		fi; # if not wallsresult[ii][1] then.. else
		#
	od;
	#
	#
	if checklevel>0 then
		check_AdjChamsdata(tweyl, torbits, wallsresult, adjdata);
	fi;
	#
	return(adjdata);
end;

