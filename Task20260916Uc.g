#Read("Task20260916Uc.g");

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



6TanCons:=function(Gram)
  local kk, trec, h, hdual, deg, tr, tinvol;
  kk:=Length(Gram)-1;
  trec:=rec(
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
  h:=MakeVectei(kk+1, 1);
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

6TanConRecs:=[];

counter:=0;
for kk in [2..10] do 
  tgras:=WeightedCompleteGraphClasses(kk);
  for tgra in tgras do 
    counter:=counter+1;
    tGram:=WeightedCompleteGraphToGram(kk, tgra);
    trec:=6TanCons(tGram);
    if trec.kk<>kk then buzz(615221); fi;
    Add(6TanConRecs, trec);
    #Printn("_______ kk=", kk, "counter", counter);
    #if not trec.nondeg then Printn("degenerate"); continue; fi;
    #if not trec.K3flag then Printn("not geom", SignatureQ(tGram)); continue; fi;
    #if not trec.isample then Printn("not ample", trec.ADEtype); continue; fi;
    if not trec.nondeg then  continue; fi;
    if not trec.K3flag then  continue; fi;
    if not trec.isample then  continue; fi;
    if not trec.3tans=[] then Printn("3tans", Length(trec.3tans)); fi;
    if Length(trec.6tans)>30 then 
      if Length(trec.6tans)>kk then 
        Printn("extra 6tans", kk, "->", Length(trec.6tans), tgra, counter); 
      fi;
    fi;
    if counter mod 100 =0 then 
      savedata(6TanConRecs);
    fi;
  od;
  savedata(6TanConRecs);
od;

#####