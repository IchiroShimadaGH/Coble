#Read("Task20260913Ua.g");


FindRatsOnK3_240:=function(GramSX, ample, uptodeg)
  local rats, extradeg, ratsdual, task, tdeg;
  #
  extradeg:=0;
  rats:=[];
  ratsdual:=[];
  #
  task:=function(tv)
    local tvdual;
    for tvdual in ratsdual do
      if tvdual*tv<0 then return(); fi;
    od;
    Add(rats, tv);
    Add(ratsdual, tv*GramSX);
    return();
  end;
  #
  tdeg:=0;
  #
  while tdeg<uptodeg do
    tdeg:=tdeg+1;
    AffES(GramSX, ample, tdeg, -2, true, task);
    if Length(rats)>10000 then return("large"); fi;
    #Printn("tdeg", tdeg, Length(rats));
  od;
  #
  return(rats);
end;


RandomXsprime:=function(kk)
  #
  local chosen, gens, newGram, newaa, newinvol, newrats, 
  trec, newrats2, tr, newtr, new6tancons;
  #
  chosen:=RandomChooseFromL(6tancons, kk);
  gens:=CopyAdd(Union(chosen), aa);
  gens:=HermiteNormalFormIntegerMat(gens);
  gens:=RemoveZeroVectsIntMat(gens);
  newGram:=TMTTmult(gens, GramS);
  newaa:=SolutionIntMat(gens, aa);
  newinvol:=RatInvolRecP2(newGram, newaa, newaa).invol;;
  if gens*invol2<>newinvol*gens then buzz(231811); fi;
  newrats:=FindRatsOnK3_240(newGram, newaa,  25);
  new6tancons:=[];
  for tr in 240rats do
    newtr:=SolutionIntMat(gens, tr);
    if newtr<>fail then Add(new6tancons, newtr); fi;
  od;
  #if not IsEqualSet(newrats2, newrats) then 
  #  buzz(44132);
  #fi;
  #
  trec:=rec(
      ran:=Length(gens),
      basis:=gens,
      newaa:=newaa,
      newGram:=newGram, 
      chosen:=chosen,
      newrats:=newrats,
      newrats:=newrats,
      newinvol:=newinvol,
      new6tancons:=new6tancons
  );
  #
  return(trec);
end;

tdatas:=[];

for tt in [1..10000] do
  for kk in [1..6] do
    trec:=RandomXsprime(kk);
    if trec.newrats="large"  then 
      Printn(tt, kk, Length(trec.newGram), "large", Length(trec.new6tancons));
    else 
      tdata:=[kk, Length(trec.newGram), Length(trec.newrats), Length(trec.new6tancons)];
      Add(tdatas, tdata);
      Printn(tt, tdata);
    fi;
  od;
  Printn("__________________");
  Printn(Collected(tdatas));
   Printn("__________________");
od;


#####