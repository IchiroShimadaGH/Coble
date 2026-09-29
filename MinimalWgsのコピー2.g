

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
    Add(EkkRecs, MakeEkkRec(Length(EkkRecs)+1));
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

SwQuasiMinimals:=function(kk, abswg)
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
  tsign, trec, tsing, flag, aa, tells, beep, givenspcons, rho, thdual;
  #
  beep:=function(beepnumb)
    localbeep("MakeK3Rec", beepnumb); Error();
  end;
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
  rho:=Length(tGram);
  givenspcons:=List([1..kk], ii->MakeVectei(kk+1, ii+1)*tGrec.phi);
  thdual:=th*tGram;
  if Set(givenspcons*thdual)<>[2] then beep(5717461); fi;
  if TMTTmult(givenspcons, tGram)<>tadj then beep(817761); fi;
  #
  trec:=rec(
    kk:=kk, 
    wg:=wg,
    adj:=tadj,
    phi:=tGrec.phi,
    Gram:=tGram,
    givenspcons:=givenspcons, 
    h:=th, 
    sign:=tsign,
    sing:=tsing, 
    ells:=tells, 
    flag:=flag
  );
  return(trec);
  #
end;



IsGkMinimal:=function(kk, wg, absstab)
  #
  local tSwgens, minflag, Swtask, ttau, newwg;
  #
  tSwgens:=GetEkkRec(kk).Swgens;
  #
  minflag:=true;
  Swtask:=function(twg, jj)
    local newtwg, tpos;
    if jj>kk then return(); fi;
    #
    Swtask(twg, jj+1);
    if minflag=false then return(); fi;
    #
    newtwg:=List(twg);
    for tpos in tSwgens[jj] do 
      newtwg[tpos]:=-newtwg[tpos];
    od;
    if newtwg<wg then minflag:=false; return(); fi;
    if newtwg<>twg then 
      Swtask(newtwg, jj+1);
      if minflag=false then return(); fi;
    fi;
    #
  end;
  #
  for ttau in absstab do 
    newwg:=ActionTauWg(ttau, wg);
    if newwg<wg then minflag:=false; break; fi;
    Swtask(List(newwg), 2);
    if not minflag then break; fi;
  od;
  return(minflag);
end;



###############

HperpRec:=function(GramS, h)
  local hdual, GramP;
  orthorec:=OrthogonalCompRec(GramS, h);
  T:=StackMats([h], orthorec.basis);
  Tinv:=InverseMat()
  return(trec);
end;

##### positive majorant

GetGramP:=function(GramS, h)
  local hdual, GramP;
  hdual:=h*GramS;
  GramP:=TransposedMat([hdual])*[hdual]-GramS;
  if h*GramP*h<>2 then beep(998912); fi;
  return(GramP);
end;

GetIniData:=function(tK3rec, GramP, basisrec)
  #
  local GramS, h, n,  ogrec, gs, hs, pos, th, tgen, thtgen, inirec;
  #
  GramS:=tK3rec.Gram;
  h:=tK3rec.h;
  n:=Length(h);
  ogrec:=OGLatFromBasisRec(basisrec);
  gs:=[IdentityMat(n)];
  hs:=[h];
  pos:=0;
  for th in hs do
    pos:=pos+1;
    for tgen in ogrec.totalgens do 
      thtgen:=th*tgen;
      if not thtgen in hs then 
        Add(hs, thtgen);
        Add(gs, gs[pos]*tgen);
      fi;
    od;
  od;
  #
  if hs<>[h, -h] then beep(44671746); fi; #Thanks to Codex
  #
  inirec:=rec(
    rank:=Length(h),
    GramS:=GramS,
    GramP:=GramP,
    basisrec:=basisrec, 
    ogrec:=ogrec,
    h:=h,
    horbit:=hs,
    transporters:=gs,
    wgs:=[tK3rec.wg],
    spconss:=[tK3rec.givenspcons]
  );
  return(inirec);
end;

IsIsomShs:=function(Sh1rec, Sh2, GramP2, basisrec2)
  #
  local beep, GramS2, h2, T, Tinv, h2Tinv,
  pos, tg;
  #
  beep:=function(beepnumb)
    localbeep("IsIsomShs", beepnumb); Error();
  end;
  #
  GramS2:=Sh2[1];
  if Sh1rec.rank<>Length(GramS2) then return(false); fi;
  h2:=Sh2[2];
  T:=IsIsomBasisRecs(Sh1rec.basisrec, basisrec2);
  if T=false then return(false); fi;
  Tinv:=InverseMat(T);
  h2Tinv:=h2*Tinv;
  if not h2Tinv in Sh1rec.horbit then 
    return(false);
  fi; 
  pos:=SinglePosition(Sh1rec.horbit, h2Tinv);
  tg:=Sh1rec.transporters[pos]*T;
  if TMTTmult(tg, GramS2)<> Sh1rec.GramS then beep(578517865); fi;
  if Sh1rec.h*tg<>h2 then beep(8121123); fi;
  return(tg);
end;


FoundShRecs:=[];

IsNewSh:=function(tK3rec)
  #
  local beep, isnewflag, nowSh, oldrec, tg, tginv, newspcons, GramP2, basisrec2;
  #
  beep:=function(beepnumb)
    localbeep("IsNewSh", beepnumb); Error();
  end;
  #
  isnewflag:=true;
  nowSh:=[tK3rec.Gram, tK3rec.h];
  GramP2:=GetGramP(tK3rec.Gram, tK3rec.h);
  basisrec2:=BasisRec(GramP2); 
  for oldrec in FoundShRecs do 
    tg:=IsIsomShs(oldrec, nowSh, GramP2, basisrec2);
    if tg<>false then 
      isnewflag:=false;
      #
      tginv:=InverseMat(tg);
      newspcons:=(tK3rec.givenspcons)*tginv;
      if TMTTmult(newspcons, oldrec.GramS)<>tK3rec.adj then beep(5817718); fi;
      Add(oldrec.spconss, newspcons);
      Add(oldrec.wgs, tK3rec.wg);
      #
      return();
    fi;
  od;
  #
  if isnewflag then 
    Add(FoundShRecs, GetIniData(tK3rec, GramP2, basisrec2));
  fi;
  #
  return();
  #
end;

######

FoundShRecs:=[];

EnhancedMinimalWgs:=function(oldminwgs, newkk)
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
        partss:=SwQuasiMinimals(newkk, twg);
       # Printn("_____", twg, "_____SwQuasiMinimals", 
       #   Product(List(partss, Length)));
        K3recs:=[];
        iterpartss:=IteratorOfCartesianProduct(partss);
        for sswg in iterpartss do 
          swg:=Flat(sswg);
          K3rec:=MakeK3Rec(newkk, swg);# MakeK3Rec is defined above.
          if K3rec.flag then 
            if IsGkMinimal(newkk, swg, stab) then 
              IsNewSh(K3rec);
              Add(K3recs, K3rec);
            fi;
          fi;
        od;
        #
        totalK3recsnops:=totalK3recsnops+ Length(K3recs);
        if newkk<8 then
          Printn(newkk, "___K3recs", Length(K3recs), totalK3recsnops, ":",
          counter, "in", nopsoldminwgs, "FoundShRecs", Length(FoundShRecs));
          Add(newminwgrecs, rec(wg:=twg, stabgens:=stabgens));
          # For newkk=8, the size of newminwgrecs would be too large.
        elif newkk<>8 then
          beep(999999999999);
        fi;
      fi;
    od; #for tav in iter do
    #
    if newkk=8 then #monitor only for newkk=8
      if counter mod 10000=0 then 
        Printn(newkk,  counter, "in", nopsoldminwgs, "FoundShRecs", Length(FoundShRecs));
      fi;
    fi;
    #
  od;#for oldminwg in oldminwgs do 
  #
  Printn("____finish", counter, "in", Length(oldminwgs), ":", 
  Length(newminwgrecs), totalK3recsnops, "FoundShRecs", Length(FoundShRecs));
  if total<>(Length(absW))^EkkRec.leng then beep(919191); fi;
  return(newminwgrecs);
 end;

NewEnhancedMinimalWgs:=function(oldminwgs, newkk) #monitor だけが上と違う
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
        partss:=SwQuasiMinimals(newkk, twg);
       # Printn("_____", twg, "_____SwQuasiMinimals", 
       #   Product(List(partss, Length)));
        K3recs:=[];
        iterpartss:=IteratorOfCartesianProduct(partss);
        for sswg in iterpartss do 
          swg:=Flat(sswg);
          K3rec:=MakeK3Rec(newkk, swg);# MakeK3Rec is defined above.
          if K3rec.flag then 
            if IsGkMinimal(newkk, swg, stab) then 
              IsNewSh(K3rec);
              Add(K3recs, K3rec);
            fi;
          fi;
        od;
        #
        totalK3recsnops:=totalK3recsnops+ Length(K3recs);
        if newkk<8 then 
          Add(newminwgrecs, rec(wg:=twg, stabgens:=stabgens));
        fi;
      fi;
    od; #for tav in iter do
    #
    if newkk<8 then
        Printn(newkk,  ":",
        counter, "in", nopsoldminwgs, "newminwgrecs", Length(newminwgrecs),  
        "totalK3recsnops", totalK3recsnops,
        "FoundShRecs", Length(FoundShRecs));
        # For newkk=8, the size of newminwgrecs would be too large.
    elif newkk=8 then #monitor only for newkk=8
      if counter mod 10000=0 then 
        Printn(newkk,  counter, "in", nopsoldminwgs, "FoundShRecs", Length(FoundShRecs));
      fi;
    else beep(587587);
    fi;
    #
  od;#for oldminwg in oldminwgs do 
  #
  Printn("____finish", counter, "in", Length(oldminwgs), ":", 
  Length(newminwgrecs), totalK3recsnops, "FoundShRecs", Length(FoundShRecs));
  if total<>(Length(absW))^EkkRec.leng then beep(919191); fi;
  return(newminwgrecs);
 end;

























#################