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
  singisotwords, oldwdbasiss;
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
  oldwdbasiss:=[[inidet, iniorec.wdbasis]];
  #
  isneworec:=function(ttdet, twdbasis)
    local torec, T, tdata;
    for tdata in oldwdbasiss do 
      if tdata[1]=ttdet  and twdbasis=tdata[2] then return(false); fi;
    od;
    return(true);
  end;
  #
  thetask:=function(torec, poss) 
    #
    # tGram=TMTTmult(tbasis, Gram); 
    # tbasisdual:=tbasis*Gram;
    # h=th*tbasis;
    #
    local pos,  tGram, th, newposs,    tw, tvdual, newtbasisdual,twdbasisinv, newwdbasis, newwdbasisinv, singflag, stw, tpos, 
    xx, newGram, ttdet, ttbasisdualinv, newaddwords, neworec, newth;
    #
    tGram:=torec.Gram;
    th:=torec.h;
    #
    for pos in poss do 
      tw:=isotwords[pos];
      newwdbasis:=AddVectToBasis(torec.wdbasis, tw);
      ttdet:=AbsInt(DeterminantIntMat(newwdbasis));
      newwdbasisinv:=InverseMat(newwdbasis);  
      singflag:=false;
      for stw in singisotwords do 
        if IsIntVect(stw*newwdbasisinv) then 
          singflag:=true;
          break; # from for stw in singisotwords do 
        fi;
      od;
      if not singflag then  
        if isneworec(ttdet, newwdbasis) then 
          Add(oldwdbasiss, [ttdet, newwdbasis]);       
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
  end;
  #
  thetask(iniorec, [1..Length(isotwords)]);
  #
  return(nonsingoverlats);
end;

