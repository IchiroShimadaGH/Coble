
#Read("Task02660927Ua.g");

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
  if Size(Group(Symkkgens))<>Factorial(kk) then buzz(651121); fi;
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

MinimalAbsWgss:=[[]];

# MinimalAbsWgs:=function(kk)
#   while Length(MinimalAbsWgss)<kk do 
#     Add(EkksRecs, MakeMinimalAbsWgs(Length(Ekks)+1));
#   od;
#   return(Ekks[kk]);
# end;
