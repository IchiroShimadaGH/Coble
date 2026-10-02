

#Read("Task20261002Pa.g");

#Random Odd Enhanced 2

Read("SplitConsTools.g");
Read("NonSingOverLats.g");
Read("ReidemeisterSchreier.g");
Read("AutShByGramP.g");


resetMT(sessionnumb);


NewIsGeom:=function(tGram, th)
  local nn;
  nn:=Length(tGram);
  if Rank(tGram)<>nn then return([false, "rank"]); fi;
  if SignatureQ(tGram)<>[nn, 1, nn-1] then return([false, "sign"]); fi;
  if th*tGram*th<>2 then return([false, "deg<>2"]); fi;
  if AffESstd(tGram, th, 0, -2, true)<>[] then return([false, "sing"]); fi;
  if AffESstd(tGram, th, 1, 0, true)<>[] then return([false, "pol"]); fi;
  if not PrimitivelyEmbeddableInK3Lattice(tGram) then return([false, "genus"]); fi;
  return(true);
end;


specialshs:=[];
savename:=Concatenation("random2Pro", String(sessionnumb));

TT:=[[-2]];

taskcounter:=0;

task:=function(adj)
  local kk, tA, ttA, newadj, trec,tGram, th,isgeom, sprats,
  cc, mm, nn, tAutShgens, tOLrecs, tOLrec, ttflag, ttGram, tth, geomflag, tAutShrec,
  geomflag0;
  #
  taskcounter:=taskcounter+1;
  #
  kk:=Length(adj);
  mm:=Maximum(50, 10*kk);
  if kk>15 then mm:=2^(kk-15)*mm;fi;
  for cc in [1..mm] do
    tA:=[RandomVectFromL(kk, [1, 3])];
    ttA:=TransposedMat(tA);
    newadj:=MatMatToMat([[adj, ttA], [tA, TT]]);
    trec:=WGraphToGramh(newadj, kk+1);
    tGram:=trec.Gram;
    th:=trec.h;
    nn:=Length(th);
    geomflag0:=NewIsGeom(tGram, th);
    if not (geomflag0=true or geomflag0=[false, "genus"]) then 
      continue;
    fi;
    #
    tAutShrec:=AutShByGramP(tGram, th);
    tAutShgens:=tAutShrec.AutShgens;
    tOLrecs:=NonSingOverLats(tGram, th, tAutShgens);
    ttflag:=false;
    for tOLrec in tOLrecs do
      ttGram:=tOLrec.Gram;
      tth:=tOLrec.h;
      geomflag:=NewIsGeom(ttGram, tth);
      if geomflag=true or geomflag=[false, "genus"] then 
         ttflag:=true;
      fi;
      if geomflag=true then 
        sprats:=GetSpRats(ttGram, tth);
        if sprats.nopss[2]> 1000 or nn>17 then 
          Add(specialshs, [sprats, [ttGram, tth]]);
          Printn(kk+1, nn, sprats.nopss, tOLrec.extdeg, ":", tAutShrec.AutShsize, ":", taskcounter);
          savedataas(specialshs, savename);
        fi;
      fi;
    od;
    if ttflag then 
      task(newadj);
    fi;
  od;
end;

task([[-2]]);



############