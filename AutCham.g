#Read("AutCham.g");

AutCham:=function(Borrec, nSaSvSss)
  local beep, vss, xx, yy, nopss, vs, Srank, basis, nowrank,
  tv, basisinv, bigvs, tvs, inibiis, ttg, bigmat, iisss,
  iiss, iis, gettg, tgss, gtask, vsset, periodaut, nonperiodaut,
  periodorbs, nonperiodorbs, dones, torb, autchamrec, idrankS;
  #
  beep:=function(beepnumb)
    localbeep("AutCham", beepnumb); Error();
  end;
  #
  vss:=List(nSaSvSss, xx->List(xx.vSliftss, yy->yy[1]));
  nopss:=List(vss, Length);
  SortParallel(nopss, vss);
  vs:=Concatenation(vss);
  #
  Srank:=Length(vs[1]);
  if Rank(vs)<>Srank then beep(55542); fi;
  basis:=[]; nowrank:=0;
  for tv in vs do 
    Add(basis, tv);
    if Rank(basis)>nowrank then 
      nowrank:=nowrank+1;
      if nowrank=Srank then 
        break; #for tv in vs do
      fi;
    else 
      if Remove(basis)<>tv then beep(43121); fi;
    fi;
  od;
  #
  basisinv:=InverseMat(basis);
  bigvs:=[];
  for tvs in vss do 
    if Intersection(tvs, basis)<>[] then 
      Append(bigvs, tvs);
    fi;
  od;
  #
  inibiis:=List(basis, tv->SinglePosition(bigvs, tv));
  bigmat:=TMTTmult(bigvs, Borrec.GramS);
  iisss:=GeneralStabChain(inibiis, bigmat);
  #
  gettg:=function(iis)
    local tg, ii;
    tg:=basisinv*List(iis, ii->bigvs[ii]);
    if TMTTmult(tg, Borrec.GramS)<>Borrec.GramS then beep(27191); fi;
    return(tg);
  end;
  #
  tgss:=[];
  for iiss in iisss do Add(tgss, List(iiss, gettg)); od;
  #
  periodaut:=[];
  nonperiodaut:=[];
  #
  vsset:=Set(vs);
  #
  idrankS:=IdentityMat(Srank);
  gtask:=function(tg, lev)
    local vstgset, tqg;
    if lev=Srank then 
      if not IsIntMat(tg) then return(true); fi;
      vstgset:=Set(vs*tg);
      if vsset<>vstgset then return(true); fi;
      Add(nonperiodaut, List(tg));
      tqg:=OLtoOqL(tg, Borrec.discS);
      if  tqg in Borrec.priodcond then 
        Add(periodaut, List(tg));
      fi;
    else 
      for ttg in tgss[lev+1] do 
        gtask(ttg*tg, lev+1); # NOT tg*tgg BUT ttg*tg
      od;
    fi;
    return(true);
  end;
  #
  gtask(idrankS, 0);
  #
  periodorbs:=[];
  dones:=[];
  for tv in vsset do 
    if not tv in dones then 
      torb:=Set(periodaut, tg->tv*tg);
      Add(periodorbs, torb);
      Append(dones, torb);
    fi;
  od;
  #
  nonperiodorbs:=[];
  dones:=[];
  for tv in vsset do 
    if not tv in dones then 
      torb:=Set(nonperiodaut, tg->tv*tg);
      Add(nonperiodorbs, torb);
      Append(dones, torb);
    fi;
  od;
  #
  autchamrec:=rec(
    #tgss:=tgss, 
    #iisss:=iisss,
    nonperiodorbs:=nonperiodorbs,
    periodorbs:=periodorbs, 
    nonperiodaut:=nonperiodaut, 
    periodaut:=periodaut
  );
  #
  return(autchamrec);
  #
end;
