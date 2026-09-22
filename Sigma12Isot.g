#Read("Sigma12Isot.g");

readdata("Sigma12rec");

GramS12:=Sigma12rec.Gram;


if DeterminantIntMat(GramS12)<>1 then beep(615261); fi;
svs:=ShortestVectors(GramS12, 2);;

if Collected(svs.norms)<>[ [ 2, 132 ] ] then buzz(51612); fi;


Read("SplitConsTools.g");

th:=MakeVectei(13, 1);
tGram:=DiagonalMats([ [[2]], (-2)*GramS12]);
tGramdual:=InverseMat(tGram);

Printn(IsGeom(tGram, th));




discrecS12:=DiscriminantForm(tGram);

nn := 13;
FFvects:=Cartesian(List([1..nn], xx->[0,1]));


Df:=2*discrecS12.discf;

MaskToVector := function(m)
    local ii, v;
    if m>=2^nn then buzz(77162); fi;
    v:=List([0..nn-1],ii -> (QuoInt(m, 2^ii) mod 2));
    return (v);
end;

VectorToMask := function(v)
  local ii, m;
  m := Sum([1..nn], ii -> v[ii] * 2^(ii-1));
  return m;
end;

MaskXor := function(a, b)
  local ans, power, i;
  ans := 0;
  power := 1;
  for i in [1..nn] do
    if (a mod 2) <> (b mod 2) then
        ans := ans + power;
    fi;
    a := QuoInt(a, 2);
    b := QuoInt(b, 2);
    power := 2 * power;
  od;
  return ans;
end;



AppendVsToBasis:=function(tbasis, vs)
  local newbasis, tv, xx;
  newbasis:=ShallowCopy(tbasis);
  Append(newbasis, vs);
  newbasis:=HermiteNormalFormIntegerMat(newbasis);
  newbasis:=Filtered(newbasis, tv-> not ForAll(tv, xx->xx=0));
  return(newbasis);
end;

FFIsotVectorsS12:=[];

for tv in FFvects do 
  if (tv*Df*tv) mod 4=0 then Add(FFIsotVectorsS12, tv); fi;
od;

FFIsotVectorsS12:=Set(FFIsotVectorsS12);

Printn("FFIsotVectorsS12", Length(FFIsotVectorsS12));
foundoverlatsS12:=[];

thetask:=function(isotspacerec)
  #
  local tv, mtv, tvrepdual,newbasisdual, newbasis, newGram,newbasisinv, newth, sings,
  newmasks, isgeom, newisotspacerec, isoflag, twdual, counter;
  #
  counter:=0;
  for tv in FFIsotVectorsS12 do
    counter:=counter+1; 
    mtv:=VectorToMask(tv);
    if mtv in isotspacerec.masks then continue; fi;
    isoflag:=true;
    for twdual in isotspacerec.addwordsdual do 
      if tv*twdual mod 2 <>0 then 
        isoflag:=false;
        break;  
      fi;
    od;
    if not isoflag then continue; fi;
    tvrepdual:=tv*discrecS12.reps_dual;
    if not IsIntVect(isotspacerec.basis*tvrepdual) then beep(66112); fi;
    newbasisdual:=AppendVsToBasis(isotspacerec.basisdual, [tvrepdual]);
    newbasis:=newbasisdual*tGramdual;
    newGram:=TMTTmult(newbasis, tGram);
    if not IsEvenLattice(newGram) then beep(813421); fi;
    newbasisinv:=InverseMat(newbasis);
    newth:=th*newbasisinv;
    if newth*newGram*newth<>2 then beep(887122); fi;
    sings:=AffESstd(newGram, newth, 0, -2, true);
    if sings<>[] then 
      #Printn("sing", isotspacerec.overdim+1, GetRootType(newGram, sings));
      continue; 
    fi;
    newmasks:=Union(isotspacerec.masks, List(isotspacerec.masks, xx->MaskXor(xx, mtv)));
    isgeom:=IsGeom(newGram, newth);
    newisotspacerec:=rec(
      masks:=newmasks, 
      basis:=newbasis, 
      basisdual:=newbasisdual, 
      Gram:=newGram, 
      th:=newth,
      isgeom:=isgeom,
      overdim:=isotspacerec.overdim+1,
      addwordsdual:=CopyAdd(isotspacerec.addwordsdual, tv*Df)
    );
    if isgeom=true then 
      Add(foundoverlatsS12, newisotspacerec);
      savedata(foundoverlatsS12);
      Printn("isgeom is true", newisotspacerec.overdim);
    fi;
    Printn("___ next level", counter, newisotspacerec.overdim, isgeom);
    thetask(newisotspacerec);
  od;
end;

iniisotspacerec:=rec(
  masks:=[0], 
  basis:=IdentityMat(nn), 
  basisdual:=tGram, 
  Gram:=tGram, 
  th:=th,
  isgeom:=IsGeom(tGram, th),
  overdim:=0,
  addwordsdual:=[]
);

thetask(iniisotspacerec);

########