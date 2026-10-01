
#Read("Task20261001Ua.g");

# Read("NonSingOverLats.g");
# Read("even_lattice_genus.g");
# Read("SplitConsTools.g");

# readdata("FoundShRecs0930Ua");

Printn("FoundShRecs0930Ua", Length(FoundShRecs0930Ua));

FoundShRecsOverLats:=[];

SpecialShRecsOverLats:=[];


GetTGenus:=function(GramS)
  local signT, signS, discT, trec;
  signS:=SignatureQ(GramS);
  if signS[2]<>1 then beep(581851); fi;
  if signS[3]=0 then beep(52231); fi;
  signT:=[22, 3, 19]-SignatureQ(GramS);
  discT:=DiscriminantForm(-GramS);
  trec:=rec(
    sign:=[signT[2], signT[3]],
    discg:= discT.discg,
    discf:= discT.discf
  );
  return(trec);
end;


counter:=0;

for trec in FoundShRecs0930Ua do
  counter:=counter+1; 
  Printn("_____________", counter);
  OLrecs:=NonSingOverLats(trec.GramS, trec.h, trec.AutShgens);
  trec.OLrecs:=OLrecs;
  Printn("__OLrecs", Length(OLrecs), List(OLrecs, xx->xx.extdeg));
  for OLrec in OLrecs do
    tGram:=OLrec.Gram;
    detS:=AbsInt(DeterminantIntMat(tGram));
    trho:=Length(tGram);
    th:=OLrec.h;
    if th*tGram*th<>2 then buzz(41764); fi;
    if AffESstd(tGram, th, 0, -2, true)<>[] then  buzz(41764); fi;
    if AffESstd(tGram, th, 1, 0, true)=[] then  
      ampleflag:=true; 
    else 
      ampleflag:=false; 
    fi;
    Tgenus:=GetTGenus(tGram);
    Trec:=EvenLatticeGenus(Tgenus.sign, Tgenus.discg, Tgenus.discf);
    #
    if Trec.count>0  then
      for tGramT in Trec.grams do 
        if SignatureQ(tGramT)<>[22-trho, 2, 20-trho] then buzz(616782); fi;
        detT:=AbsInt(DeterminantIntMat(tGramT));
        if detT<>detS then buzz(571121); fi;
      od;
    fi;
    #
    sprats:=GetSpRats(tGram, th);
    sconsinttypes:=SpconsIntTypes(tGram, sprats.scons);
    sprats.sconstypes:=sconsinttypes;
    Printn("_____",  "rho", trho, ampleflag, ": Tcount", Trec.count, 
    ": sprats", sprats.nopss, sconsinttypes);
    OLrec.rho:=trho;
    OLrec.ampleflag:=ampleflag;
    OLrec.Trec:=Trec;
    OLrec.sprats:=sprats;
    OLrec.sprats:=sprats;
    #
    if ampleflag and Trec.count>0 then 
      Add(FoundShRecsOverLats, OLrec);
    else 
      Add(SpecialShRecsOverLats, OLrec);
    fi;
  od;
  if counter mod 20=0 then 
    savedataas(FoundShRecsOverLats, "tempFoundShRecsOverLats"); 
    savedataas(SpecialShRecsOverLats, "tempSpecialShRecsOverLats");
  fi;
  Printn(counter, Length(FoundShRecsOverLats), Length(SpecialShRecsOverLats));
od;

savedata(FoundShRecsOverLats);
savedata(SpecialShRecsOverLats);


########## 

