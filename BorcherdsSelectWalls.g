#Read("Walls.g");

IsWallBestCoeffs:=function(vs, w)
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
	pivotfunc, NB;
	#
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
			if printlevel>1 then
				Printn("pivot check done: this should be removed; counter=", counter);
			fi;
		fi;
		#
	od;
	#
	Printn("too many pivots");
	return(IsWallMinSubscript(vs, w));
	#
end;


check_Walls:=function(orbits, result)
	#
	#
	local beep, walls, nonwalls, wallsamples,
	nonwallsamples, ii, res,
	tw, pos, newres, x0, c, vs, v, coeffs, ttw;
	#
	beep:=function(beepn) Printn("wrong check_Walls", beepn); Error(); end;
	#
	if Length(result)<>Length(orbits) then
		beep(0);
	fi;
	#
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
			newres:=IsWallBestCoeffs(walls, tw);
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
	if printlevel>1 then
		Printn("check_Walls done: this should be removed.");
	fi;
	return(true);
end;

Make_Walls:=function(orbits)
	local result, ii, torbit, w, vs, jj, st, res, checktrial;
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
			res:=IsWallBestCoeffs(vs, w);
		fi;
		Add(result, res);
	od;
	if checklevel>0 then
		check_Walls(orbits, result);
	fi;
	return(result);
end;


check_AdjChamsdata:=function(oldweyl, torbits, wallsresult, adjdata)
	#
	local beep, ii, tr, newweyl;
	#
	beep:=function(beepn) Printn("wrong in check_AdjChamsdata", beepn); Error(); end;
	if Length(adjdata)<>Length(torbits) then
		beep(0);
	fi;
	for ii in [1..Length(adjdata)] do
		if adjdata[ii][1]=_nw then
			if wallsresult[ii][1]<>false then beep(1); fi;
		elif adjdata[ii][1]=_mtc then
			tr:=adjdata[ii][2];
			if tr*GramSX*tr<>-2 then beep(2); fi;
			if tr*GramSX*ampleX<=0 then beep(3); fi;
		elif adjdata[ii][1]=_newcham then
			newweyl:=adjdata[ii][2];
			#
			if checklevel>3 and not IsWeylVect(GramL, newweyl) then 
				beep(4);
			fi;
		fi;
	od;
	if printlevel>1 then
		Printn("check_AdjChamsdata done: this should be removed;");
	fi;
	return(true);
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
		Printn("wrong BorcherdsSelectWalls", numb); Error();
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
