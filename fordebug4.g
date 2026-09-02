
GeneralStabChain:=function(inibiis, bigmat)
  #
  local beep, nn,NN,  targetmat, ii, jj, orbits, orbit,
  candidatess, candidates,  extendtask, reachflag, foundpiis, 
  aii, bii, lev, inipiis, 
  xx, yy, kk, tv;
  #
  beep:=function(beepnumb)
    localbeep("GeneralStabChain", beepnumb);Error();
  end;
  #
  nn:=Length(inibiis);
  NN:=Length(bigmat);
  if bigmat<>TransposedMat(bigmat) then beep(68215873); fi;
  targetmat:=List(inibiis, ii->List(inibiis, jj->bigmat[ii][jj]));
  candidatess:=[];
  for jj in [1..nn] do
    bii:=inibiis[jj];
    aii:=bigmat[bii][bii];
    candidates:=[];
    for kk in [1..NN] do 
      if bigmat[kk][kk]=aii then
        Add(candidates, kk); 
      fi;
    od;
    Add(candidatess, candidates);
  od;
  #
  reachflag:=false;
  foundpiis:=[];
  extendtask:=function(piis, leng)
    local targetv, yy, yyv, fflag, pos, tii;
    if leng=nn then 
      reachflag:=true; 
      foundpiis:=List(piis);
    else 
      targetv:=targetmat[leng+1];
      for yy in candidatess[leng+1] do 
        yyv:=bigmat[yy];
        fflag:=true;
        pos:=0;
        for tii in piis do 
          pos:=pos+1;
          if yyv[tii]<>targetv[pos] then 
            fflag:=false; 
            break; #from for tii in piis
          fi;
        od;
        if fflag then 
          Add(piis, yy);
          extendtask(piis, leng+1);
          if Remove(piis)<>yy then beep(999121); fi;
          if reachflag then 
            break; # from for yy in candidatess[leng+1] do 
          fi;
        fi;
      od;
    fi;
    return(true);
  end;
  #
  orbits:=[];
  for lev in [1..nn] do
    inipiis:=List([1..lev-1], jj->inibiis[jj]);
    orbit:=[];
    for xx in candidatess[lev] do 
      reachflag:=false;
      Add(inipiis, xx);
      extendtask(inipiis, lev);
      if Remove(inipiis)<>xx then beep(9919); fi;
      if reachflag then 
        Add(orbit, foundpiis);
      fi;
      
    od;
    Add(orbits, orbit);
  od;
  #
  return(orbits);
end;

AutDiscfByGeneralStabChain:=function(arg...)
  local discv, discg, discf, vs, bigmat, orbits, gens, gensperm, orbit,
  piis, tg, kk, ii;
  if Length(arg)=2 then 
    discg:=arg[1];
    discf:=arg[2];
  elif Length(arg)=1 then 
    discg:=arg.discg;
    discf:=arg.discf;
  else beep(33911);
  fi;
  vs:=Cartesian(List(discg, kk->[0..kk-1]);
  bigmat:=CopyNormalDiscf(TMTTmult(vs, discf));
  orbits:=GeneralStabChain(discf, bigmat);
  #
  gens:=[];
  gensperm:=[];
  #
  discv:=function(vv)
    local dvv, pos, xx;
    dvv:=[];
    pos:=0;
    for xx in vv do 
      pos:=pos+1;
      Add(dvv, xx mod discg[pos]);
    od;
    return(dvv);
  end;
  for orbit in orbits do
    for piis in orbit do 
      tg:=List(piis, ii->vs[ii]);
      if discf<>CopyNormalDiscf(TMTTmult(tg, discf)) then 
        beep(787711)
      fi;
    od; 
    Add(gens, tg);
    tgp:=PermList(List(List(vs*tg, discv), tw->Position(vs, tw)));
    Add(gensperm, tgp);
  od;
  #
  size:=
  #
end;
