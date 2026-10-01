#Read("NonSingOverLats.g");

AddVectToBasis:=function(basis, vv)
  local newbasis, xx;
  newbasis:=ShallowCopy(basis); 
  Add(newbasis, vv);
  newbasis:=HermiteNormalFormIntegerMat(newbasis);
  newbasis:=Filtered(newbasis, xx->not IsZeroVect(xx));
  return(newbasis);
end;


NonSingOverLats:=function(Gram, h, autShgens)
  #
  local nn, Gramdual, discrec, discf, thdual, discg,
  iter, isotwordsdual, tw, twdual, iniorec, nonsingoverlats, 
  isneworec, thetask, isotwords, overlats, inidet, extdeg,
  tvdual, newtbasisdual, newth, newGram, wdbasis,poss, 
  singisotwords, oldwdbasiss, lengdiscg, moddiscg, diagmatdiscg, 
  autShgensperms, tg, tqg,isotwordstg, theG,
  IsGminimal,  GetMinimalInGorb, GetAllElements, GetOLrec,
  nonsingOLrecs, minelmss;
  #
  nn:=Length(h);
  Gramdual:=InverseMat(Gram);
  discrec:=DiscriminantForm(Gram);
  discf:=discrec.discf;
  discg:=discrec.discg;
  lengdiscg:=Length(discg);
  moddiscg:=function(wd)
    return(List([1..lengdiscg], ii-> (wd[ii] mod discg[ii])));
  end;
  diagmatdiscg:=DiagonalMat(discg);
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
  autShgensperms:=[];
  for tg in autShgens do 
    tqg:=OLtoOqL(tg, discrec);
    isotwordstg:=List(isotwords, tw->moddiscg(tw*tqg));
    poss:=List(isotwordstg, tw->Position(isotwords, tw));
    Add(autShgensperms, PermList(poss));
  od;
  #
  theG:=Group(autShgensperms);
  Printn("__size of G", Size(theG));
  #
  IsGminimal:=function(iis)
    local tgp, iistgp;
    for tgp in theG do 
      iistgp:=List(iis, ii->ii^tgp);
      Sort(iistgp);
      if iistgp<iis then return(false); fi;
    od;
    return(true);
  end;
  #
  GetMinimalInGorb:=function(iis)
    local tgp, iistgqp, min;
    min:=iis;
    for tgp in theG do 
      iistgqp:=List(iis, ii->ii^tgp);
      Sort(iistgqp);
      if iistgqp<min then 
        min:=iistgqp;
      fi;
    od;
    return(min);
  end;
  #
  GetAllElements:=function(wdbasisinv)
    local iis, ii, tw;
    iis:=[];
    ii:=0;
    for tw in isotwords do 
      ii:=ii+1;
      if IsIntVect(tw*wdbasisinv) then 
        Add(iis, ii);
      fi;
    od;
    return(iis);
  end;
  #
  inidet:=AbsInt(DeterminantIntMat(Gram));
  #
  GetOLrec:=function(minelmiis)
    #
    local extgroupgens, newbasisdual, newth, newGram, nopsxtgroup, OLrec;
    #
    extgroupgens:=CopyAppend(diagmatdiscg, List(minelmiis, ii->isotwords[ii]));
    extgroupgens:=HermiteNormalFormIntegerMat(extgroupgens);
    extgroupgens:=Filtered(extgroupgens, tw->not IsZeroVect(tw));
    #
    newbasisdual:=CopyAppend(Gram, extgroupgens*discrec.reps_dual);
    newbasisdual:=HermiteNormalFormIntegerMat(newbasisdual);
    newbasisdual:=Filtered(newbasisdual, tw->not IsZeroVect(tw));
    #
    newth:=SolutionIntMat(newbasisdual, thdual);
    newGram:=TMTTmult(newbasisdual, Gramdual);
    #
    nopsxtgroup:=Length(minelmiis)+1; # +1 is for the zero word.
    if nopsxtgroup^2*AbsInt(DeterminantIntMat(newGram))<>inidet then beep(919111); fi;
    #
    OLrec:=rec(
      Gram:=newGram, 
      h:=newth, 
      basisdual:=newbasisdual,
      extgroupgens:=extgroupgens, 
      extdeg:=nopsxtgroup
    );
    return(OLrec);
  end;
  #
  iniorec:=rec(
    addwords:=[],
    addwordiis:=[],
    wdbasis:=List(diagmatdiscg), 
    wdbasisinv:=InverseMat(diagmatdiscg)
  );
  #
  minelmss:=[[]]; # the zero wd is not included in isotwords 
  nonsingOLrecs:=[GetOLrec([])];
  #
  #
  thetask:=function(torec, poss) 
    #
    #
    local pos, newposs, oldiis, newiis, 
     tw, tvdual, newtbasisdual,twdbasisinv, newwdbasis, 
     newwdbasisinv, singflag, stw, tpos, hashnewwdbasis, singflag2, ddinv, 
    newOLrec, newelms, minnewelms, 
    xx, newGram, ttdet, ttbasisdualinv, newaddwords, neworec, newth, dd; ;
    #
    oldiis:=torec.addwordiis;
    #
    for pos in poss do 
      newiis:=CopyAdd(oldiis, pos);
      if not IsGminimal(newiis) then continue; fi; #in for pos in poss do 
      #
      tw:=isotwords[pos];
      newwdbasis:=AddVectToBasis(torec.wdbasis, tw);
      newwdbasisinv:=InverseMat(newwdbasis);  
      singflag:=false;
      for stw in singisotwords do 
        if IsIntVect(stw*newwdbasisinv) then 
          singflag:=true;
          break; # from for stw in singisotwords do 
        fi;
      od;
      if singflag then continue; fi; #in for pos in poss do 
      #
      newelms:=GetAllElements(newwdbasisinv);
      minnewelms:=GetMinimalInGorb(newelms);
      if not minnewelms in minelmss then 
        Add(minelmss, minnewelms);
        newOLrec:=GetOLrec(minnewelms);
        Add(nonsingOLrecs, newOLrec);
        Printn("__new overlattice", newOLrec.extdeg, Length(nonsingOLrecs));
      fi;
      #
      newaddwords:=ShallowCopy(torec.addwords); Add(newaddwords, tw);
      neworec:=rec(
        wdbasis:=newwdbasis,
        wdbasisinv:=newwdbasisinv,
        addwords:=newaddwords,
        addwordiis:=newiis
      );
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
    od; #for pos in poss do 
    #
  end;
  #
  thetask(iniorec, [1..Length(isotwords)]);
  #
  return(nonsingOLrecs);
end;

