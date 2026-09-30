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
  singisotwords, oldwdbasiss, hash, maxgg, lev1counter;
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
  #
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
  #
  maxgg:=Maximum(discg);
  hash:=function(wdmat)
    local fmat, aa, xx;
    fmat:=Flat(wdmat);
    aa:=0;
    for xx in fmat do 
      aa:=aa*maxgg+xx;
    od;
    return(aa);
  end;
  #
  oldwdbasiss:=Set([hash(iniorec.wdbasis)]);
  #
  lev1counter:=0;
  #
  thetask:=function(torec, poss) 
    #
    # tGram=TMTTmult(tbasis, Gram); 
    # tbasisdual:=tbasis*Gram;
    # h=th*tbasis;
    #
    local pos,  tGram, th, newposs,    tw, tvdual, newtbasisdual,twdbasisinv, newwdbasis, newwdbasisinv, singflag, stw, tpos, hashnewwdbasis, singflag2, ddinv, 
    xx, newGram, ttdet, ttbasisdualinv, newaddwords, neworec, newth, dd; #, st, Rtm1, Rtm2;
    #
    tGram:=torec.Gram;
    th:=torec.h;
    #
    for pos in poss do 
      tw:=isotwords[pos];
      newwdbasis:=AddVectToBasis(torec.wdbasis, tw);
      ttdet:=AbsInt(DeterminantIntMat(newwdbasis));
      newwdbasisinv:=InverseMat(newwdbasis);  
      #st:=Runtime();
      dd:=Lcm(List(Flat(newwdbasisinv), DenominatorRat));
      ddinv:=dd*newwdbasisinv;
      singflag:=false;
      for stw in singisotwords do 
        if IsZeroVect((stw*ddinv) mod dd) then 
          singflag:=true;
          break; # from for stw in singisotwords do 
        fi;
      od;
      # Rtm1:=Runtime()-st;
      # st:=Runtime();
      # singflag2:=false;
      # for stw in singisotwords do 
      #   if IsIntVect(stw*newwdbasisinv) then 
      #     singflag2:=true;
      #     break; # from for stw in singisotwords do 
      #   fi;
      # od;
      # Rtm2:=Runtime()-st;
      # Printn(Rtm1, Rtm2);
      # if singflag<>singflag2 then beep(66152); fi;
      if not singflag then  
        hashnewwdbasis:=hash(newwdbasis);
        if not hashnewwdbasis  in oldwdbasiss then 
          AddSet(oldwdbasiss, hashnewwdbasis);       
          tvdual:=tw*discrec.reps_dual;
          newtbasisdual:=AddVectToBasis(torec.basisdual, tvdual);
          newth:=SolutionIntMat(newtbasisdual, thdual);
          newGram:=TMTTmult(newtbasisdual, Gramdual);
          newaddwords:=ShallowCopy(torec.addwords); Add(newaddwords, tw);
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
          Add(nonsingoverlats, neworec);
          #Printn("___", Length(nonsingoverlats), Length(newaddwords), 
          #Length(oldwdbasiss));
          newposs:=[];
          for tpos in poss do 
            if tpos>pos then 
              if IsInt(tw*isotwordsdual[tpos]) then 
                if not IsIntVect(isotwords[tpos]*newwdbasisinv) then 
                  Add(newposs, tpos);
                fi;
              fi;
            fi;
          od;
          thetask(neworec, newposs);
        fi;
      fi;
    od; 
    #
    if Length(torec.addwords)=1 then 
      lev1counter:=lev1counter+1;
      Printn("___in NonSingOverLats__", lev1counter, "in", Length(isotwords), 
            ":", Length(nonsingoverlats));
    fi;
    #
  end;
  #
  thetask(iniorec, [1..Length(isotwords)]);
  #
  return(nonsingoverlats);
end;

