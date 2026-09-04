#Read("Adj.g");


IsWeylVector:=function(Borrec, tweyl) # used in check_AdjChamsdata
	local mG, GramL;
  GramL:=Borrec.GramL;
	if not IsPrimitiveVect(tweyl) then return(false); fi;
	if not tweyl*GramL*tweyl=0 then return(false); fi;
	mG:=(-1)*(HoloSphereLattice(GramL, tweyl).Gram);
	if Length(ShortestVectors(mG, 2).vectors)>0 then return(false); fi;
	return(true);
end;

IsWalldefv:=function(Borrec, defv) 
  local defvdual;
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
  local beep, GramL, nrmv, mm, rLs, aa, defvdual, 
  aadefvL, aanrmv, nRvRsrec,
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
    oldweyl:=List(tweyl),
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
  IsWeylVector(Borrec, adjrec.newweyl);
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
  if orthoroots<>Union(rLs, (-1)*rLs) then beep(42211); fi;
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
    parats:=List(ttrLs, xr->(u*GramL*xr)/(adjrec.oldweyl*GramL*xr));
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

###################