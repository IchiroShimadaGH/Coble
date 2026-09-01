# Read("uvabdiscf.g");



_uumat:=[[0,1/2], [1/2,0]];
_vvmat:=[[1,1/2], [1/2, 1]];
_aamat:=[[1/2]];
_bbmat:=[[3/2]];



isisomuvab:=function(uvab1, uvab2)
  local uvab1s, pos, tuvab, tp, tpp, dduvab, nnuvab, xx, 
  iseven1, iseven2, leng1, leng2;
  if uvab1=uvab2 then return(true); fi;
  #
  leng1:=[2,2,1,1]*uvab1;
  leng2:=[2,2,1,1]*uvab2;
  #
  if leng1<>leng2 then return(false); fi;
  #
  iseven1:=(uvab1[3]=0) and (uvab1[4]=0);
  iseven2:=(uvab2[3]=0) and (uvab2[4]=0);
  #
  if iseven1<>iseven2 then return(false); fi;
  #
  #
  uvab1s:=[uvab1];
  pos:=1;
  while pos<=Length(uvab1s) do 
    tuvab:=uvab1s[pos];
    for tp in  uvab_isompairs do
      tpp:=tp;
      dduvab:=tuvab-tpp[1];
      if Set(dduvab, xx->xx>=0)=[true] then 
        nnuvab:=dduvab+tpp[2];
        if nnuvab=uvab2 then return(true);fi;
        if not nnuvab in uvab1s then Add(uvab1s, nnuvab); fi;
      fi;
       tpp:=Reversed(tp);
       dduvab:=tuvab-tpp[1];
       if Set(dduvab, xx->xx>=0)=[true] then 
        nnuvab:=dduvab+tpp[2];
        if nnuvab=uvab2 then return(true);fi;
        if not nnuvab in uvab1s then Add(uvab1s, nnuvab); fi;
      fi;
    od;
    pos:=pos+1;
  od;
  return(false);
end;




uvabvect:=function(discrec)
  local ntqfrec, ntqf, xx, tuvab;
  if Set(discrec.discg)<>[2] then beep(52875872); fi;
  ntqfrec:=Normalform_fqf(discrec.discg, discrec.discf);
  ntqf:=ntqfrec.decompdata[1].pfqf;
  tuvab:=List([_uumat, _vvmat, _aamat, _bbmat], xx->OccurNumb(ntqf, xx));
  return(tuvab);
end;

uvab2qf:=function(tuvabvect)
  local dmatslist, uvabmatslist, ii, kk;
  dmatslist:=[];
  uvabmatslist:=[_uumat, _vvmat, _aamat, _bbmat];
  ii:=0;
  for kk in tuvabvect do
    ii:=ii+1;
    Append(dmatslist, List([1..kk], xx->uvabmatslist[ii]));
  od;
  return(DiagonalMats(dmatslist));
end;

#AutDiscfByStabilizerChain:=function(inibiis, bigmat)

uvabOqLOrder:=function(tuvabvect)
  local tqf, leng, vs, jj, inibiis, tv, bigmat, tsize;
  tqf:=uvab2qf(tuvabvect);
  leng:=[2,2,1,1]*tuvabvect;
  vs:=Cartesian(List([1..leng], jj->[0,1]));
  inibiis:=List(IdentityMat(leng), tv->SinglePosition(vs, tv));
  bigmat:=CopyNormalDiscf(TMTTmult(vs, tqf));
  tsize:=AutDiscfByStabilizerChain(inibiis, bigmat).order;
  return(tsize);
end; 







