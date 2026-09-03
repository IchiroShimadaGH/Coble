#Read("Walls.g");

#####################
# 2026/09/0: For Coble.
#
# There are a collection of complex codes (main: IsWallBestCoeffs and  check_Walls) and 
# a simple code (BorcherdsSelectWalls).
#
###################

IsWallMinSubscript:=function(Borrec, vs, w)
	#
	# this will be used in IsWallBestCoeffs,
	# when IsWallBestCoeffs gives up after too many pivots.
	#
	# We assume that vs is the full rank.
	# Does the hyperplane w*x=0 intersect with the interior of
	# the cone defined by v*x>=0 for all v in vs?
	# This is a problem of linear programming.
	# The answer is yes if and only if the problem
	#       minimize w*x
	#       subject to v*x>=0
	# has unbounded to -infinity,
	# or if and only if there exists x0 such that
	# w*x0<0 and v*x0>=0 for any v in vs,
	# or if and only if w can NOT  be written as
	# a linear combination of some v's in v with
	# positive coefficients.
	# If true, then it returns [true, x0].
	# If false, then it returns [false, [coeffs, v's]].
	#
	local dic, f, n, nonbasics, svsiis, NBinv, #NBinv is nonbasic inv;
	counter, entcand, pos, c, tflag, y,evpos, sol,
	entvarnumb, entpos, ynumb, leavevarnumb, leavey,a,correcty, pivotfunc,
	NB, rankSX;
	#
  rankSX:=Length(Borrec.GramS);
	n:=Length(w);
	dic:=CopyMat(vs);
	f:=CopyVect(w);
	#
	# initialize
	#
	svsiis:=FindSubvsWithRank(vs, rankSX);
	if svsiis=false then Printn("cannot apply IsWall", RankMat(vs)); Error(); fi;
	nonbasics:=svsiis[2];
	NBinv:=Inverse(svsiis[1]);
	Apply(dic, v->v*NBinv);
	f:=f*NBinv;
	#
	counter:=0;
	#
	while counter<10*Length(vs) do
		#
		counter:=counter+1;
		#
		# find entering
		#
		entcand:=[];
		pos:=0; # position in nonbasiscs
		for c in f do
			pos:=pos+1;
			if c<0 then
				#
				tflag:=true;
				for y in dic do
					if y[pos]<0 then tflag:=false; break; fi;
				od;
				if tflag then #unbound!
					NBinv:=Inverse(List(nonbasics, xx->vs[xx]));
					evpos:=MakeVectei(rankSX, pos);
					sol:=NBinv*evpos;
					#Printn(counter, "unbound");
					return([true, sol]);
				fi;
				#
				Add(entcand, nonbasics[pos]);
			fi;
		od;
		#
		if Length(entcand)=0 then # bound
			#Printn(counter, "bound");
			return([false, [f, List(nonbasics, ii->vs[ii])]]);
		fi;
		#
		entvarnumb:=Minimum(entcand);
		entpos:=Position(nonbasics, entvarnumb);
		#
		# find leaving
		#
		ynumb:=0;
		for y in dic do
			ynumb:=ynumb+1;
			if y[entpos]<0 then leavevarnumb:=ynumb; leavey:=y; break; fi;
		od;
		#
		# pivot
		#
		nonbasics[entpos]:=leavevarnumb;
		a:=leavey[entpos];
		correcty:=CopyVect(leavey);
		correcty[entpos]:=-1;
		correcty:=(-1/a)*correcty;
		pivotfunc:=function(yy)
			local aa;
			aa:=yy[entpos];
			yy[entpos]:=0;
			yy:=yy+aa*correcty;
			return(yy);
		end;
		Apply(dic, pivotfunc);
		f:=pivotfunc(f);
		#
		#  check
		if counter <0 then
			NB:=List(nonbasics, ii->vs[ii]);
			if dic*NB<>vs then Printn("wrong pivot 1"); Error(); fi;
			if f*NB<>w then Printn("wrong pivot 2"); Error(); fi;
			#if printlevel>1 then
			#	Printn("pivot check done: this should be removed; counter=", counter);
			#fi;
		fi;
    #
	od;
	#
	Printn("too many pivots"); Error();
	#
end;


IsWallBestCoeffs:=function(Borrec, vs, w)
	#
	# We assume that vs is the full rank and that
	# the region defined by w*x>=0 and v*x>=0 for all v in vs
	# has an interior point.
	# Does the hyperplane w*x=0 intersect with the interior of
	# the cone defined by v*x>=0 for all v in vs?
	# This is a problem of linear programming.
	# The answer is yes if and only if the problem
	#       minimize w*x
	#       subject to v*x>=0
	# has unbounded to -infinity,
	# or if and only if there exists x0 such that
	# w*x0<0 and v*x0>=0 for any v in vs,
	# or if and only if w can NOT  be written as
	# a linear combination of some v's in vs with
	# positive coefficients.
	# If true, then it returns [true, x0].
	# If false, then it returns [false, [coeffs, v's]].
	#
	local dic, f, n, nonbasics, svsiis, NBinv,
	#NBinv is nonbasic inv;
	counter, minc,  y,evpos, sol,
	entvarnumb, entpos, ynumb, leavevarnumb, leavey,a,correcty,
	pivotfunc, NB, rankSX;
	#
  rankSX:=Length(Borrec.GramS);
	n:=Length(w);
	dic:=CopyMat(vs);
	f:=CopyVect(w);
	#
	# initialize
	#
	svsiis:=FindSubvsWithRank(vs, rankSX);
	if svsiis=false then
		Printn("cannot apply IsWall", RankMat(vs)); Error();
	fi;
	nonbasics:=svsiis[2];
	NBinv:=Inverse(svsiis[1]);
	Apply(dic, v->v*NBinv);
	f:=f*NBinv;
	#
	counter:=0;
	#
	while counter<1000 do
		#
		counter:=counter+1;
		#
		# find entering
		#
		minc:=Minimum(f);
		#
		if minc>=0 then # bound
			return([false, [f, List(nonbasics, ii->vs[ii])]]);
		fi;
		#
		entpos:=Position(f, minc);
		entvarnumb:=nonbasics[entpos];
		#
		# find leaving
		#
		ynumb:=0;
		a:=0;
		for y in dic do
			ynumb:=ynumb+1;
			if y[entpos]<a then
				a:=y[entpos];
				leavevarnumb:=ynumb;
				leavey:=y;
			fi;
		od;
		#
		if a=0 then
			NBinv:=Inverse(List(nonbasics, xx->vs[xx]));
			evpos:=MakeVectei(rankSX, entpos);
			sol:=NBinv*evpos;
			return([true, sol]); #Printn(counter, "unbound");
		fi;
		#
		# pivot
		#
		nonbasics[entpos]:=leavevarnumb;
		correcty:=CopyVect(leavey);
		correcty[entpos]:=-1;
		correcty:=(-1/a)*correcty;
		pivotfunc:=function(yy)
			local newy, aacorrecty;
			aacorrecty:=yy[entpos]*correcty;
			yy[entpos]:=0;
			newy:=yy+aacorrecty;
			return(newy);
		end;
		Apply(dic, pivotfunc);
		f:=pivotfunc(f);
		#
		#  check
		if counter <0 then
			NB:=List(nonbasics, ii->vs[ii]);
			if dic*NB<>vs then Printn("wrong pivot 1"); Error(); fi;
			if f*NB<>w then Printn("wrong pivot 2"); Error(); fi;
			#if printlevel>1 then
			#	Printn("pivot check done: this should be removed; counter=", counter);
			#fi;
		fi;
		#
	od;
	#
	Printn("too many pivots");
	return(IsWallMinSubscript(Borrec, vs, w));
	#
end;


check_Walls:=function(Borrec, orbits, result)
	#
	#
	local beep, walls, nonwalls, wallsamples,
	nonwallsamples, ii, res,
	tw, pos, newres, x0, c, vs, v, coeffs, ttw, rankSX;
	#
	beep:=function(beepn) localbeep("wrong check_Walls", beepn); Error(); end;
	#
	if Length(result)<>Length(orbits) then
		beep(0);
	fi;
	#
  rankSX:=Length(Borrec.GramS);
	walls:=[];
	nonwalls:=[];
	wallsamples:=[];
	nonwallsamples:=[];
	ii:=0;
	#
	for res in result do
		ii:=ii+1;
		if res[1] then
			Append(walls, orbits[ii]);
			Add(wallsamples, RandomOp(orbits[ii]));
		else
			Append(nonwalls, orbits[ii]);
			Add(nonwallsamples, RandomOp(orbits[ii]));
		fi;
	od;
	if RankMat(walls)<>rankSX then beep(1); fi;
	#
	for tw in wallsamples do
		pos:=Position(walls, tw);
		Remove(walls, pos);
		#
		# tw is NOT in walls from now
		#
		if RankMat(walls)=rankSX then
			# Note that, if RankMat(walls)<rankSX, then tw is a wall.
			newres:=IsWallBestCoeffs(Borrec, walls, tw);
			if not newres[1] then beep(2); fi;
			x0:=newres[2];
			if tw*x0>=0 then beep(3); fi;
			for ttw in walls do
				if ttw*x0<0 then beep(4); fi;
			od;
		fi;
		#
		# tw returns to walls
		#
		Add(walls, tw);
	od;
	#
	#
	for tw in nonwallsamples do
		newres:=IsWallBestCoeffs(walls, tw);
		if newres[1] then beep(5); fi;
		coeffs:=newres[2][1];
		vs:=newres[2][2];
		if coeffs*vs<>tw then beep(6); fi;
		for c in coeffs do if c<0 then beep(7); fi; od;
		for v in vs do if not v in walls then beep(8); fi;od;
	od;
	#
	#if printlevel>1 then
	Printn("__check_Walls done: this should be removed.");
	#fi;
	return(true);
end;

Make_Walls:=function(Borrec, orbits)
	local result, ii, torbit, w, vs, jj, st, res, checktrial, rankSX;
  #
  rankSX:=Length(Borrec.GramS);
	result:=[];
	ii:=0;
	for torbit in orbits do
		ii:=ii+1;
		w:=torbit[1];
		vs:=[];
		for jj in [1..Length(orbits)] do
			if ii<jj or (jj<ii and result[jj][1]) then
				Append(vs, orbits[jj]);
			fi;
		od;
		if RankMat(vs)<rankSX then
			res:=[true, "rank reason"];
		else
			res:=IsWallBestCoeffs(Borrec, vs, w);
		fi;
		Add(result, res);
	od;
	#if checklevel>0 then
	#	check_Walls(orbits, result);
	#fi;
  check_Walls(Borrec, orbits, result);
	return(result);
end;




#########################



BorcherdsSelectWalls:=function(tlinforms)
	#
	# The input should be normalized as
	# tlinforms=Set(tlinforms, lf->GetPrimitiveIntVect(lf));
	# and does not contain a zerovector.
	#
	# The list of elements of tlinforms
	# that define a wall of the cone defined by x*tfl>=0 for all tfl in tlinforms.
	#
	local beep, dim, tlf, walllfs, lf, wallflag,  st, ntlinforms, sst2, ccc,
	cPrintn, n2tlinforms;
	#
	beep:=function(numb)
		localbeep("wrong BorcherdsSelectWalls", numb); Error();
	end;
	#
	#
	#if globalprintlevel>0 then
  if -1>0 then
		cPrintn:=function(arg) Printn(arg); end;
	else
		cPrintn:=function(arg)  end;
	fi;
	#
	ntlinforms:=Set(tlinforms, GetPrimitiveIntVect);
	n2tlinforms:=Set(tlinforms, GetPrimitiveIntVect);
	#
	dim:=Length(ntlinforms[1]);
	#
	sst2:=Runtime();
	#cPrintn("_____start BorcherdsSelectWalls", Length(tlinforms));
	#
	if Length(ntlinforms)=dim then
		#
		# tlinforms defines the cone over a simplex. #Added 2018 01/21
		#
		#cPrintn("_____done BorcherdsSelectWalls", Length(walllfs), StringTime(Runtime()-sst2));
		return(ntlinforms);
	fi;
	#
	walllfs:=[];
	ccc:=0;
	for tlf in n2tlinforms do
		#
		# Since we use SubtractSet and AddSet to ntlinforms,
		# we cannot write "for tlf in ntlinforms do" #Added 2018/07/21
		#
		ccc:=ccc+1;
		SubtractSet(ntlinforms, [tlf]); #Add tfl back below.
		if Rank(ntlinforms)<dim then
			#
			# if removed then the rest is not bounded;
			# hence tlf defines a wall. #Added 2018 01/21
			#
			wallflag:=[true];
		else
			#IsAWall:=function(vs, w, checklevel) #in LatticeTools
			wallflag:=IsAWall(ntlinforms, tlf, 1);
		fi;
		if wallflag[1] then
			Add(walllfs, tlf);
			AddSet(ntlinforms, tlf); #ntlinforms should be a SET.
			#
			# No need to put tlf back if not wallflag[1]
			#
		fi;
		#
		if ccc mod 1000 =0 then
			cPrintn("_____in BorcherdsSelectWalls", ccc, "in", Length(tlinforms), Length(walllfs), StringTime(Runtime()-sst2));
		fi;
	od;
	#
	#cPrintn("_____done BorcherdsSelectWalls", Length(walllfs), StringTime(Runtime()-sst2));
	#
	walllfs:=Set(walllfs);
	return(walllfs);
end;
