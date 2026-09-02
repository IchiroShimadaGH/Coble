ttAutDiscfByGeneralStabChain:=function(arg)
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
    Printn("disc is trivial");
    autqrec:=rec(
      discg:=discg,
      discf:=discf,
      gens:=[], 
      gensperm:=[],
      size:=1
    );
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
  #
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
  #
  autqrec:=rec(
    discg:=discg,
    discf:=discf,
    gens:=gens, 
    gensperm:=gensperm,
    size:=size,
    inibiis:=inibiis,
    bigmat:=bigmat, 
    orbits:=orbits
  );
  #
  return(autqrec);
  #
end;

