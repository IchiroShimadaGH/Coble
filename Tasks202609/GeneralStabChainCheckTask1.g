#Read("GeneralStabChainCheckTask1.g");

Read("GeneralStabChain.g");


for tt in [1..10000] do
  for kk in [1..2] do  
    L:=Union([-kk..kk], [0,0,0,0,0,0,0,0]);
    for dim in [2..6] do 
      tGram:=2*RandomEvenLatticeFromL(dim, L);
      tdet:=AbsInt(DeterminantIntMat(tGram));
      if tdet>3000 or tdet=1 then continue; fi;# in  for dim in [2..10] do 
      tdisc:=DiscriminantForm(tGram);
      if tdisc.discg=0 then continue; fi;
      tautqrec:=AutDiscfByGeneralStabChain(tdisc);
      for uu1 in [1..2] do 
        tU:=RandomUnimodMat(dim);
        ttGram:=TMTTmult(tU, tGram);
        ttdisc:=DiscriminantForm(ttGram);
        ttautqrec:=AutDiscfByGeneralStabChain(ttdisc);
        if tautqrec.size<>ttautqrec.size then beep(776645); fi;
      od;
      for uu2 in [1..2] do 
        tU:=RandomUnimodMat(dim);
        ttGram:=TMTTmult(tU, tGram);
        ttdisc:=DiscriminantForm(ttGram);
        ttautqrec:=AutDiscfByGeneralStabChain(ttdisc.discg, ttdisc.discf);
        if tautqrec.size<>ttautqrec.size then beep(736645); fi;
      od;
      Printn(tt, "kk", kk, "dim", dim, "disc", 
      tdisc.discg, "size", tautqrec.size);
    od;
  od;
od;



