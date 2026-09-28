Read("IsNewSh.g");

#positive majorant

GetGramP:=function(GramS, h)
  local hdual, GramP;
  hdual:=h*GramS;
  GramP:=TransposedMat([hdual])*[hdual]-GramS;
  if h*GramP*h<>2 then beep(998912); fi;
  return(GramP);
end;

GetIniData:=function(GramS, h)
  #
  local n, GramP, ogrec, gs, hs, pos, th, tgen, thtgen, inirec;
  #
  n:=Length(h);
  GramP:=GetGramP(GramS);
  basisrec:=BasisRec(GramP); 
  ogrec:=OGLatFromBasisRec(basisrec);
  gs:=[IdentityMat(n)];
  hs:=[h];
  pos:=0;
  for th in hs do
    pos:=pos+1;
    for tgen in ogrec.totalgens do 
      thtgen:=th*tgen;
      if not thtgen in hs then 
        Add(hs, thtgen);
        Add(gs, gs[pos]*tgen);
      fi;
    od;
  od;
  inirec:=rec(
    rank:=Length(h),
    GramS:=GramS,
    GramP:=GramP,
    basisrec:=basisrec, 
    ogrec:=ogrec,
    h:=h,
    horbit:=hs,
    transporters:=gs
  );
  return(inirec);
end;

IsIsomSh:=function(Sh1rec, Sh2pair)
  GramS2:=Sh2pair[1];
  if Sh1rec.rank<>Length( GramS2) then return(false); fi;
  h2:=Sh2pair[2];
  GramP2:=GetGramP(GramS2, h2);
  basisrec2:=BasisRec(GramP2); 
end;

######