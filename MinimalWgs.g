

#Read("MinimalWgs.g");

Read("SplitConsTools.g");

MakeEkk:=function(kk)
  local Ekk, ii, jj;
  Ekk:=[];
  for jj in [2..kk] do 
    for ii in [1..jj-1] do 
      Add(Ekk, [ii, jj]);
    od;
  od;
  return(Ekk);
end;


Ekkpos:=function(kk, iijj)
  local ii, jj;
  ii:=iijj[1]; jj:=iijj[2];
  if ii>=jj then beep(665522); fi;
  return((jj-1)*(jj-2)/2+ii);
end;

Ekkiijj:=function(kk, pos)
  local ii, jj, aa;
  aa:=1+8*pos;
  jj:=2 + QuoInt(RootInt(8*pos - 7) - 1, 2);
  ii:=pos-(jj-1)*(jj-2)/2;
  if ii>=jj then beep(72919); fi;
  return([ii, jj]);
end;

TauToLargetau:=function(kk, Ekk, tau)
  local iijj, xx, poss, largetau;
  poss:=List(Ekk, iijj->Ekkpos(kk, CopySort(List(iijj, xx->xx^tau))));
  largetau:=PermList(poss);
  return(largetau);
end;

MakeEkkRec:=function(kk)
  local trec, Symkkgens, aa, tau, poss, largetau, Swkkgens, xx, Ekk;
  #
  Ekk:=MakeEkk(kk);
  #
  Symkkgens:=[];
  for aa in [1..kk-1] do 
    largetau:=TauToLargetau(kk, Ekk, (aa, aa+1));
    Add(Symkkgens, largetau);
  od;
  if Size(Group(Symkkgens))<>Factorial(kk) then 
    if kk<>2 then 
      buzz(651121); 
    fi;
  fi;
  #
  Swkkgens:=[];
  for xx in [1..kk] do 
    Add(Swkkgens, Filtered([1..Length(Ekk)], tpos->xx in Ekk[tpos]));
  od;
  #
  trec:=rec(
    kk:=kk, 
    Ekk:=Ekk,
    leng:=Length(Ekk),
    Symkkgens:=Symkkgens, 
    Swgens:=Swkkgens
  );
  return(trec);
end;

EkkRecs:=[[]];

GetEkkRec:=function(kk)
  while Length(EkkRecs)<kk do 
    Add(EkkRecs, MakeEkkRec(kk));
  od;
  return(EkkRecs[kk]);
end;




ActionTauWg:= function(largetau, wg)
  local newposs, newwg, tpos;
  newposs:=List([1..Length(wg)], tpos->tpos^largetau);
  newwg:=List(newposs, tpos->wg[tpos]);
  return(newwg);
end;

absW:=[0,1,2];

MinimalAbsWgss:=[[]];


GetGeneratingSetSmall:=function(totalelms)
  local targetsize, newgensset, tsize, tg, ttsize;
  targetsize:=Length(totalelms);
  newgensset:=[];
  tsize:=1;
  while tsize<targetsize do 
    tg:=Random(totalelms);
    Add(newgensset, tg);
    ttsize:=Size(Group(newgensset));
    if ttsize=tsize then 
      if Remove(newgensset)<>tg then buzz(47147); fi;
    else 
      tsize:=ttsize;
    fi; 
  od;
  if tsize<>targetsize then buzz(47111); fi;
  return(newgensset);
end;


MakeMinimalAbsWgs:=function(oldminwgs, newkk)
  local oldkk,  EkkRec, newminwgs, oldminwg, iter, ii,
  twg, tav, minflag, tg, tgtwg, counter;
  #
  oldkk:=newkk-1;
  EkkRec:=GetEkkRec(newkk);
  newminwgs:=[];
  counter:=0;
  for oldminwg in oldminwgs do 
    counter:=counter+1;
    iter:=IteratorOfCartesianProduct(List([1..newkk-1], ii->absW));
    twg:=List(oldminwg);
    for tav in iter do
      twg:=CopyAppend(oldminwg, tav);
      minflag:=true;
      for tg in Group(EkkRec.Symkkgens) do 
        tgtwg:=ActionTauWg(tg, twg);
        if tgtwg<twg then minflag:=false; break; fi;
      od; 
      if minflag then 
        Add(newminwgs, twg);
      fi;
    od;
    if counter mod 100=0 then 
      Printn("____", counter, "in", Length(oldminwgs), ":", Length(newminwgs));
    fi;
  od;
  Printn("____finish", counter, "in", Length(oldminwgs), ":", Length(newminwgs));
  return(newminwgs);
 end;

theWs:=[-2,-1,0,1,2];

SwMinimals:=function(kk, abswg)
  local tabs, jj, tpos, partss, ttab, tintval,
  parts, subwg, nzposs,nznops, task, tabpos;
  #
  tabs:=List([1..kk], jj->jj*(jj-1)/2);
  if Length(abswg)<>tabs[kk] then buzz(812811); fi;
  tpos:=1;
  partss:=[];
  for tabpos in [1..Length(tabs)-1] do
    tintval:=[tabs[tabpos]+1..tabs[tabpos+1]];
    if tintval<>[] then 
      parts:=[];
      subwg:=List(tintval, tpos->abswg[tpos]);
      nzposs:=Filtered([1..Length(tintval)], ttpos-> subwg[ttpos]<>0);
      if nzposs=[] then 
        parts:=[List(subwg)];
      else 
        nznops:=Length(nzposs);
        subwg[nzposs[1]]:=-subwg[nzposs[1]];
        #
        parts:=[];
        task:=function(tsubwg, ff)
          local ffpos;
          if ff>nznops then 
            Add(parts,List(tsubwg));
          else 
            ffpos:=nzposs[ff];
            task(tsubwg, ff+1);
            tsubwg[ffpos]:=-tsubwg[ffpos];
            task(tsubwg, ff+1);
          fi;
        end;
        #
        task(subwg, 2);
      fi;
      Add(partss, parts);
    fi;
  od;
  return(partss);
end;

MakeK3Rec:=function(kk, wg)
  local wwg, xx, tadj, ttadj, pos, iijj,   tGrec, tGram, th, 
  tsign, trec, tsing, flag, aa, tells;
  #
  if Length(wg)<>kk*(kk-1)/2 then buzz(471471); fi;
  #
  wwg:=List(wg, xx->xx+2);
  tadj:=(-2)*IdentityMat(kk);
  pos:=0;
  for aa in wwg do 
    pos:=pos+1;
    iijj:=Ekkiijj(kk, pos);
    tadj[iijj[1]][iijj[2]]:=aa;
    tadj[iijj[2]][iijj[1]]:=aa;
  od; 
  ttadj:=Cover22(tadj);
  tGrec:=ModLatticeKerRec(ttadj);
  tGram:=tGrec.redGram;
  th:=MakeVectei(kk+1, 1)*tGrec.phi;
  if th*tGram*th<>2 then beep(18821); fi;
  tsign:=SignatureQ(tGram);
  flag:=true;
  if tsign[2]<>1 then 
    tsing:="N/A";
    tells:="N/A";
    flag:=false; 
  else  
    tsing:=AffESstd(tGram, th, 0, -2, true);
    if tsing<>[] then 
      flag:=false;  
    fi;
    tells:=AffESstd(tGram, th, 1, 0, true);
    if tells<>[] then 
      flag:=false; 
    fi;
  fi; 
  #
  trec:=rec(
    kk:=kk, 
    wg:=wg,
    adj:=tadj,
    Ladj:=ttadj,
    phi:=tGrec.phi,
    Gram:=tGram,
    h:=th, 
    sign:=tsign,
    sing:=tsing, 
    ells:=tells, 
    flag:=flag
  );
  return(trec);
  #
end;






Enhanced2MakeMinimalAbsWgs:=function(oldminwgs, newkk)
  local oldkk,  EkkRec, newminwgrecs, oldminwg, iter, ii,
  twg, tav, minflag, tg, tgtwg, counter, stab, total,
  orbsize, Gksize, stabgens, partss, iterpartss, swg, sswg,
  K3rec,  K3recs, totalK3recsnops, nopsoldminwgs;
  #
  oldkk:=newkk-1;
  EkkRec:=GetEkkRec(newkk);
  Gksize:=Size(Group(EkkRec.Symkkgens));
  newminwgrecs:=[];
  counter:=0;
  total:=0;
  totalK3recsnops:=0;
  nopsoldminwgs:=Length(oldminwgs);
  for oldminwg in oldminwgs do 
    counter:=counter+1;
    iter:=IteratorOfCartesianProduct(List([1..newkk-1], ii->absW));
    twg:=List(oldminwg);
    for tav in iter do
      twg:=CopyAppend(oldminwg, tav);
      stab:=[];
      minflag:=true;
      for tg in Group(EkkRec.Symkkgens) do 
        tgtwg:=ActionTauWg(tg, twg);
        if tgtwg<twg then minflag:=false; break; fi;
        if tgtwg=twg then Add(stab, tg); fi;
      od; 
      if minflag then 
        stabgens:=GetGeneratingSetSmall(stab);
        orbsize:=Gksize/Length(stab);
        if not IsInt(orbsize) then buzz(163671); fi;
        total:=total+orbsize;
        #
        partss:=SwMinimals(newkk, twg);
       # Printn("_____", twg, "_____SwMinimals", 
       #   Product(List(partss, Length)));
        K3recs:=[];
        iterpartss:=IteratorOfCartesianProduct(partss);
        for sswg in iterpartss do 
          swg:=Flat(sswg);
          K3rec:=MakeK3Rec(newkk, swg);
          if K3rec.flag then 
            Add(K3recs, K3rec);
          fi;
        od;
        #
        
        totalK3recsnops:=totalK3recsnops+ Length(K3recs);
        Printn(newkk, "___K3recs", Length(K3recs), totalK3recsnops, ":",
        counter, "in", nopsoldminwgs);
        Add(newminwgrecs, rec(wg:=twg, stabgens:=stabgens, K3recs:=K3recs));
      fi;
    od;
  od;
  Printn("____finish", counter, "in", Length(oldminwgs), ":", 
  Length(newminwgrecs), totalK3recsnops);
  if total<>(Length(absW))^EkkRec.leng then beep(919191); fi;
  return(newminwgrecs);
 end;



























#################