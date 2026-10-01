

#Read("Task20260924Ub.g");

Printn("wrong"); beep(581758712); 
Read("NonSingOverLats.g");

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

Ekks:=[[]];

GetEkk:=function(kk)
  while Length(Ekks)<kk do 
    Add(Ekks, MakeEkk(Length(Ekks)+1));
  od;
  return(Ekks[kk]);
end;


Getss:=function(mask, kk)
  local ss, ii;
  ss := List([1..kk], ii -> 1);
  for ii in [1..kk] do
      if QuoInt(mask, 2^(ii-1)) mod 2 = 1 then
          ss[ii] := -1;
      fi;
  od;
  return(ss);
end;




IsMinimal:=function(kk, wg)
  local stab, Ekk, tau, poss, tauwg, mask, ss, pos,
  newtauwg, aa, iijj, stabs;
  #
  stabs:=[];
  Ekk:=GetEkk(kk);
  for tau in SymmetricGroup(kk) do
    poss:=List(Ekk, iijj->Ekkpos(kk, CopySort([iijj[1]^tau, iijj[2]^tau])));
    tauwg:=List(poss, pos->wg[pos]);
    for mask in [0..2^kk-1] do 
      ss:=Getss(mask, kk);
      pos:=0;
      newtauwg:=[];
      for aa in tauwg do
        pos:=pos+1;
        iijj:=Ekkiijj(kk,pos);
        if ss[iijj[1]]*ss[iijj[2]]=-1 then 
          Add(newtauwg, 4-aa);
        else 
          Add(newtauwg, aa);
        fi;
      od;
      tauwg:=newtauwg;
      if tauwg<wg then return(false); fi;
      if tauwg=wg then Add(stabs, [mask, tau]); fi;
    od;
  od;
  #
  return([true, stabs]);
end;


Makeaddvs:=function(kk)
  local addvs, task;
  addvs:=[];
  #
  task:=function(jj, addv)
    local aa, minflag;
    if jj=kk then 
      Add(addvs, List(addv));
    else 
      for aa in [0..4] do 
        Add(addv, aa);
        minflag:=true;
        #
        if List(addv, xx->4-xx)<addv then 
          # the action of ss[jj]=1, ss[kk+1]=-1, tau=id gives a smaller element.
          minflag:=false;
        fi;
        #
        if minflag then 
          task(jj+1, addv);
        fi;
        if aa<>Remove(addv) then beep(221221); fi; 
      od;

    fi;
  end;
  #
  task(0, []);
  #
  return(addvs);
end;



NewWGraphToGramh:=function(w, kk)
  #
  local aa, iijj, ii, jj, tv, tGram, redGrec, redGram,
  th, thdual, gensplcons, M, trec, xx, pos;
  #
  if Length(w)<>kk*(kk-1)/2 then buzz(471471); fi;
  #
  M:=(-2)*IdentityMat(kk);
  pos:=0;
  for aa in w do 
    pos:=pos+1;
    iijj:=Ekkiijj(kk, pos);
    M[iijj[1]][iijj[2]]:=aa;
    M[iijj[2]][iijj[1]]:=aa;
  od; 
  #
  # the order of w is different from the original WGraphToGramh.
  #
  tGram:=Cover22(M);
  redGrec:=ModLatticeKerRec(tGram);
  redGram:=redGrec.redGram;
  th:=MakeVectei(kk+1, 1)*redGrec.phi;
  thdual:=th*redGram;
  if th*thdual<>2 then beep(18821); fi;
  #
  #
  gensplcons:=[];
  for jj in [1..kk] do 
    tv:=MakeVectei(kk+1, jj+1)*redGrec.phi;
    if tv*redGram*tv<>-2 then beep(88111); fi;
    if tv*thdual<>2 then beep(281834); fi;
    Add(gensplcons, tv);
  od;
  if TMTTmult(gensplcons, redGram)<>M then beep(916621); fi;
  #
  #
  trec:=rec(
    kk:=kk, 
    adj:=M, 
    h:=th, 
    gensplcons:=gensplcons, 
    Gram:=redGram
  );
  return(trec);
end;


uptokk:=10;
savename:="testresult";

TotalRecs:=[];

thetask:=function(kk, wg)
  local addvs, addv, newwg, trec, tGram, th, isgeom, minflag, scons,
  nooverlatflag, t0Gram, t0h, orecs, oreccounter, ttrec, sprats;
  addvs:=Makeaddvs(kk);
  for addv in addvs do 
    newwg:=CopyAppend(wg, addv);
    #koko(1);
    trec:=NewWGraphToGramh(newwg, kk+1);#Here we must use New version!
    t0Gram:=trec.Gram;
    t0sign:=SignatureQ(t0Gram);
    if t0sign[2]<>1 then continue; fi; # in for addv in addvs do 
    t0h:=trec.h;
    orecs:=NonSingOverLats(t0Gram,  t0h);
    Printn("kk", kk, "orecs", Length(orecs));
    oreccounter:=0;
    for orec in orecs do
      oreccounter:=oreccounter+1;
      ttrec:=StructuralCopy(trec);
      ttrec.orec:=orec;
      tGram:=orec.Gram;
      th:=orec.h;
      isgeom:=IsGeom(tGram, th);
      ttrec.isgeom:=isgeom;
      # ttmdiscrec:=DiscriminantForm(-tGram);
      # tsign:=[22, 3, 19]-SignatureQ(tGram);
      # ttsign:=[tsign[2], tsign[3]];
      # if tsign[3]>0 then 
      #   Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
      # fi;
      if isgeom=true then
        sprats:=GetSpRats(tGram, th);
        # if Trec.count=0 then buzz(47652);fi;
        # if Trec.count>1 then Printn("____ more than one T", Trec.count); fi;
        ttrec.sprats:=sprats;
        # trec.Trec:=Trec;
        Printn("kk", Length(newwg), "rho", Length(tGram), Collected(newwg),  
            "extdeg", orec.extdeg, "sprats", sprats.nopss,   
            ":",  oreccounter, "in", Length(orecs)); 
        #
        Add(TotalRecs, ttrec);
        savedataas(TotalRecs, savename);
        if kk+1<=uptokk then 
          thetask(kk+1, newwg);
        fi;
      else
        # if tsign[3]>0 and .count<>0 then buzz(222652);fi;
        Printn("kk", Length(newwg), "rho", Length(tGram), Collected(newwg),   
        "not isgeom", isgeom,
        ":",  oreccounter, "in", Length(orecs));
      fi;
      #
     
      
    od;
  od;
end;

uptokk:=3;
savename:="minrecs20260924";
TotalRecs:=[];

thetask(1, []);


############