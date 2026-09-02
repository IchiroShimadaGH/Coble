#Read("GeneralStabChain.g");

GeneralStabChain:=function(inibiis, bigmat)
  #
  local beep, nn,NN,  targetmat, ii, jj, orbits, orbit,
  candidatess, candidates,  extendtask, reachflag, foundpiis, 
  aii, bii, lev, inipiis, sscandidatess, sscandidates,
  xx, yy, kk, tv, vyy, targetv, fflag, fingerprints, pos;
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
  sscandidatess:=[];
  for jj in [1..nn] do
    if nn=1 then 
      sscandidates:=List(candidatess[1]);
    else 
      sscandidates:=[];
      targetv:=targetmat[jj];
      for yy in candidatess[jj] do 
        vyy:=bigmat[yy];
        if targetv[jj]<>vyy[yy]  then beep(33211); fi; #for check 
        fflag:=true;
        for pos in [1..jj-1] do 
          if targetv[pos]<>vyy[inibiis[pos]] then
            fflag:=false; break; #from for pos in [1..nn-1] do 
          fi;
        od;
        if fflag then 
          Add(sscandidates, yy);
        fi;
      od;
    fi;
    Add(sscandidatess, sscandidates);
  od;
  #
  fingerprints:=List(sscandidatess, Length);
  #
  reachflag:=false;
  foundpiis:=[];#just a holder
  extendtask:=function(piis, leng)
    local targetv, yy, yyv, fflag, pos, tii, piismat,
    newsscandidates;
    if leng=nn then 
      reachflag:=true; 
      piismat:=List(piis, ii->List(piis, jj->bigmat[ii][jj])); #for check
      if piismat<>targetmat then beep(31198); fi; #for check
      foundpiis:=List(piis);
    else 
      targetv:=targetmat[leng+1];
      newsscandidates:=[];
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
          Add(newsscandidates, yy);
        fi;
      od;
      if Length(newsscandidates)<>fingerprints[leng+1] then 
        return(true);
      fi;
      for yy in newsscandidates do 
        Add(piis, yy);
        extendtask(piis, leng+1);
        if Remove(piis)<>yy then beep(999121); fi;
        if reachflag then 
          break; # from  for yy in newsscandidates do
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
    for xx in sscandidatess[lev] do 
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

AutDiscfByGeneralStabChain:=function(arg)
  local beep, discv, discg, discf, vs, bigmat, orbits, gens, gensperm, orbit,
  piis, tg, kk, ii, size, autqrec, tw, tgp, leng, inibiis, tv;
  #
  #
  beep:=function(beepnumb)
    localbeep("AutDiscfByGeneralStabChain", beepnumb);Error();
  end;
  #
  if Length(arg)=2 then 
    discg:=arg[1];
    discf:=arg[2];
  elif Length(arg)=1 then 
    discg:=arg[1].discg;
    discf:=arg[1].discf;
  else beep(33911);
  fi;
  if discg=[] then 
    autqrec:=rec(
      discg:=discg,
      discf:=discf,
      gens:=[], 
      gensperm:=[],
      size:=1
    );
    return(autqrec);
  fi;
  vs:=Cartesian(List(discg, kk->[0..kk-1]));
  bigmat:=CopyNormalDiscf(TMTTmult(vs, discf));
  leng:=Length(discg);
  inibiis:=List(IdentityMat(leng), tv->Position(vs, tv));
  orbits:=GeneralStabChain(inibiis, bigmat);
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
        beep(787711);
      fi;
      AddSet(gens, tg);
      tgp:=PermList(List(List(vs*tg, discv), tw->Position(vs, tw)));
      if List(inibiis, xx->xx^tgp)<>piis then beep(767687); fi;
      AddSet(gensperm, tgp);
    od; 
  od;
  #
  size:=Size(Group(gensperm));
  if size<>Product(List(orbits, Length)) then beep(88811); fi;
  #
  autqrec:=rec(
    discg:=discg,
    discf:=discf,
    gens:=gens, 
    gensperm:=gensperm,
    size:=size
  );
  #
  return(autqrec);
  #
end;
