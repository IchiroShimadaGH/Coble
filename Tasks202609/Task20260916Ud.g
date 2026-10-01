#Read("Task20260916Ud.g");

Read("WeightedCompleteGraphClasses.g");

NpriFindRatsOnK3:=function(GramSX, ample, uptodeg)
  local rats, extradeg, ratsdual, task, tdeg;
  #
  extradeg:=0;
  rats:=[];
  ratsdual:=[];
  tdeg:=0;
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
  while tdeg<uptodeg do
    tdeg:=tdeg+1;
    AffES(GramSX, ample, tdeg, -2, true, task);
    #Printn("tdeg", tdeg, Length(rats));
  od;
  #
  return(rats);
end;



New6TanCons:=function(Gram, h)
  local kk, trec,  hdual, deg, tr, tinvol;
  kk:=Length(Gram)-1;
  trec:=rec(
    h2:=h,
    kk:=kk,
    nondeg:=true,
    K3flag:=false,
    isample:=false,
    discrec:="void", 
    sings:="void",
    ADEtype:=0, 
    invol:="void",
    rats2:="void",
    3tans:="void",
    6tans:="void",   
  );
  hdual:=h*Gram;
  if Rank(Gram)<>kk+1 then  trec.nondeg:=false; fi;
  if trec.nondeg then 
    if SignatureQ(Gram)=[kk+1, 1, kk] then 
      trec.K3flag:=PrimitivelyEmbeddableInK3Lattice(Gram);
    fi;
    if trec.K3flag then 
      trec.discrec:=DiscriminantForm(Gram);
      trec.sings:=AffESstd(Gram, h, 0, -2, true);
      if trec.sings=[] then 
        trec.isample:=true;
        trec.3tans:=[];
        trec.6tans:=[];
        tinvol:=RatInvolRecP2(Gram, h, h).invol;
        trec.invol:=tinvol;
        trec.rats2:=NpriFindRatsOnK3(Gram, h,  2);
        for tr in  trec.rats2 do 
          deg:=tr*hdual;
          if deg=1 then 
            if tr*tinvol=tr then buzz(51162); fi;
            Add(trec.3tans, Set([tr*tinvol, tr]));
          elif deg=2 then 
            if tr*tinvol=tr then buzz(52262); fi;
            AddSet(trec.6tans, Set([tr*tinvol, tr]));
          else 
            buzz(7162711);
          fi;
        od;
      else 
        trec.ADEtype:=GetRootType(Gram, trec.sings);
      fi;
    fi;
  fi;
  return(trec);
end;

tGram:=DiagonalMat([2,-2,-2,-2,-2,-2,-2,-2,-2]);
th:=[3,-1,-1,-1,-1,-1,-1,-1,-1];

ttrec:=New6TanCons(tGram, th);;


#####