#Read("NonSingOverLats.g");

AddVectToBasis:=function(basis, vv)
  local newbasis, xx;
  newbasis:=ShallowCopy(basis); 
  Add(newbasis, vv);
  newbasis:=HermiteNormalFormIntegerMat(newbasis);
  newbasis:=Filtered(newbasis, xx->not IsZeroVect(xx));
  return(newbasis);
end;


NonSingOverLats:=function(Gram, h)
  #
  local nn, Gramdual, discrec, discf, thdual, discg,
  iter, isotwordsdual, tw, twdual, iniorec, nonsingoverlats, 
  isneworec, thetask, isotwords, overlats, inidet, extdeg,
  tvdual, newtbasisdual, newth, newGram, wdbasis,
  singisotwords, oldwdbasisinvs;
  #
  nn:=Length(h);
  Gramdual:=InverseMat(Gram);
  discrec:=DiscriminantForm(Gram);
  discf:=discrec.discf;
  discg:=discrec.discg;
  thdual:=h*Gram;
  #
  iter:=IteratorOfCartesianProduct(List(discrec.discg, ii->[0..ii-1]));
  isotwords:=[];
  isotwordsdual:=[];
  singisotwords:=[];

  for tw in iter do 
    if IsZeroVect(tw) then continue; fi;
    twdual:=tw*discf;
    if modtZ(twdual*tw)=0 then 
      tvdual:=tw*discrec.reps_dual;
      newtbasisdual:=AddVectToBasis(Gram, tvdual); 
      newth:=SolutionIntMat(newtbasisdual, thdual);
      newGram:=TMTTmult(newtbasisdual, Gramdual);
      if AffESstd(newGram, newth, 0, -2, true)=[] then 
        Add(isotwordsdual, twdual);
        Add(isotwords, tw);
      else 
        Add(singisotwords, tw);
      fi;
    fi;
  od;
  #
  Printn("__nonsingisotwords", Length(isotwords));
  Printn("__singisotwords", Length(singisotwords));
  #
  inidet:=Product(discg);
  iniorec:=rec(
    Gram:=Gram, 
    h:=h, 
    basisdual:=Gram, 
    addwords:=[],
    det:=inidet,
    extdeg:=1,
    wdbasis:=DiagonalMat(discg), 
    wdbasisinv:=InverseMat(DiagonalMat(discg))
  );
  #
  nonsingoverlats:=[iniorec];
  oldwdbasisinvs:=[[inidet, iniorec.wdbasisinv]];
  #
  isneworec:=function(ttdet, twdbasis)
    local torec, T, tdata;
    for tdata in oldwdbasisinvs do 
      if tdata[1]=ttdet then 
        T:=twdbasis*tdata[2];
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
    local pos, twdual, tGram, th, tbasisdual, taddwords,   tw, tvdual, newtbasisdual,twdbasisinv, newwdbasis, newwdbasisinv, singflag, stw, 
    xx, newGram, ttdet, ttbasisdualinv, newaddwords, neworec, newth;
    #
    tGram:=torec.Gram;
    th:=torec.h;
    tbasisdual:=torec.basisdual;
    taddwords:=torec.addwords;
    twdbasisinv:=torec.wdbasisinv;
    #
    pos:=0;
    for twdual in isotwordsdual do 
      pos:=pos+1;
      if taddwords=[] or ForAll(taddwords*twdual, IsInt) then 
        tw:=isotwords[pos];
        if not IsIntVect(tw*twdbasisinv) then 
          newwdbasis:=AddVectToBasis(torec.wdbasis, tw);
          ttdet:=AbsInt(DeterminantIntMat(newwdbasis));
          if isneworec(ttdet, newwdbasis) then 
            newwdbasisinv:=InverseMat(newwdbasis);
            Add(oldwdbasisinvs, [ttdet, newwdbasisinv]);
            singflag:=false;
            for stw in singisotwords do 
              if IsIntVect(stw*newwdbasisinv) then 
                singflag:=true;
                break; # from for stw in singisotwords do 
              fi;
            od;
            if not singflag then 
              tvdual:=tw*discrec.reps_dual;
              newtbasisdual:=AddVectToBasis(tbasisdual, tvdual);
              newth:=SolutionIntMat(newtbasisdual, thdual);
              newGram:=TMTTmult(newtbasisdual, Gramdual);
              newaddwords:=ShallowCopy(taddwords); Add(newaddwords, tw);
              neworec:=rec(
                Gram:=newGram, 
                h:=newth, 
                basisdual:=newtbasisdual, 
                det:=ttdet,
                extdeg:=inidet/ttdet,
                wdbasis:=newwdbasis,
                wdbasisinv:=newwdbasisinv,
                addwords:=newaddwords
              );
              if neworec.extdeg<0 then beep(88811); fi;
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

