#Read("NonSingOverLats.g");

NonSingOverLats:=function(Gram, h)
  #
  local nn, Gramdual, discrec, discf, thdual, 
  iter, isotwordsdual, tw, twdual, iniorec, nonsingoverlats, 
  isneworec, thetask, isotwords, overlats, inidet, extdeg;
  #
  nn:=Length(h);
  Gramdual:=InverseMat(Gram);
  discrec:=DiscriminantForm(Gram);
  discf:=discrec.discf;
  thdual:=h*Gram;
  #
  iter:=IteratorOfCartesianProduct(List(discrec.discg, ii->[0..ii-1]));
  isotwords:=[];
  isotwordsdual:=[];

  for tw in iter do 
    if IsZeroVect(tw) then continue; fi;
    twdual:=tw*discf;
    if modtZ(twdual*tw)=0 then 
      Add(isotwordsdual, twdual);
      Add(isotwords, tw);
    fi;
  od;
  #
  Printn("__isotwords", Length(isotwords));
  #
  inidet:=DeterminantIntMat(Gram);
  iniorec:=rec(
    Gram:=Gram, 
    h:=h, 
    basisdual:=Gram, 
    addwords:=[],
    det:=inidet,
    extdeg:=1
  );
  #
  nonsingoverlats:=[iniorec];
  overlats:=[iniorec];
  #
  isneworec:=function(ttdet, ttbasisdualinv)
    local torec, T;
    for torec in overlats do 
      if torec.det=ttdet then 
        T:=torec.basisdual*ttbasisdualinv;
        if IsIntMat(T) then return(false); fi;
      fi;
    od;
    return(true);
  end;
  #
  thetask:=function(torec) 
    #
    # tGram=TMTTmult(tbasis, Gram); 
    # tbasisdual:=tbasis*Gram;
    # h=th*tbasis;
    #
    local pos, twdual, tGram, th, tbasisdual, taddwords,   tw, tvdual, newtbasisdual,
    tT, H, xx, newGram, ttdet, ttbasisdualinv, newaddwords, neworec, newth;
    #
    tGram:=torec.Gram;
    th:=torec.h;
    tbasisdual:=torec.basisdual;
    taddwords:=torec.addwords;
    #
    pos:=0;
    for twdual in isotwordsdual do 
      pos:=pos+1;
      if taddwords=[] or ForAll(taddwords*twdual, IsInt) then 
        tw:=isotwords[pos];
        tvdual:=tw*discrec.reps_dual;
        if SolutionIntMat(tbasisdual, tvdual)=fail then 
          tT:=ShallowCopy(tbasisdual); Add(tT, tvdual);
          H:=HermiteNormalFormIntegerMat(tT);
          newtbasisdual:=Filtered(H, xx->not IsZeroVect(xx));
          newth:=SolutionIntMat(newtbasisdual, thdual);
          newGram:=TMTTmult(newtbasisdual, Gramdual);
          ttdet:=DeterminantIntMat(newGram);
          ttbasisdualinv:=InverseMat(newtbasisdual);
          if isneworec(ttdet, ttbasisdualinv) then 
            neworec:=rec(
              Gram:=newGram, 
              h:=newth, 
              basisdual:=newtbasisdual, 
              det:=ttdet
            );
            Add(overlats, neworec);
            if AffESstd(newGram, newth, 0, -2, true)=[] then 
              newaddwords:=ShallowCopy(taddwords); Add(newaddwords, tw);
              neworec.addwords:=newaddwords;
              extdeg:=SqRoot(inidet/ttdet);
              neworec.extdeg:=extdeg;
              #
              if not IsInt(extdeg) then beep(777111); fi;
              Add(nonsingoverlats, neworec);
              thetask(neworec);
            fi;
          fi;
        fi;
      fi;
    od; 
    #
  end;
  #
  thetask(iniorec);
  #
  return(nonsingoverlats);
end;

