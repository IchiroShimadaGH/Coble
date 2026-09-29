#Read("IsNewSh.g");

##### positive majorant

GetGramP:=function(GramS, h)
  local hdual, GramP;
  hdual:=h*GramS;
  GramP:=TransposedMat([hdual])*[hdual]-GramS;
  if h*GramP*h<>2 then beep(998912); fi;
  return(GramP);
end;

GetIniData:=function(tK3rec)
  #
  local GramS, h, n, GramP, ogrec, gs, hs, pos, th, tgen, thtgen, inirec, basisrec;
  #
  GramS:=tK3rec.Gram;
  h:=tK3rec.h;
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
    transporters:=gs,
    wgs:=[tK3rec.wg],
    spconss:=[tK3rec.givenspcons]
  );
  return(inirec);
end;

IsIsomShs:=function(Sh1rec, Sh2)
  #
  local beep, GramS2, h2, GramP2, basisrec2, T, Tinv, h2Tinv,
  pos, tg;
  #
  beep:=function(beepnumb)
    localbeep("IsIsomShs", beepnumb); Error();
  end;
  #
  GramS2:=Sh2[1];
  if Sh1rec.rank<>Length(GramS2) then return(false); fi;
  h2:=Sh2[2];
  GramP2:=GetGramP(GramS2, h2);
  basisrec2:=BasisRec(GramP2); 
  T:=IsIsomBasisRecs(Sh1rec, basisrec2);
  if T=false then return(false); fi;
  Tinv:=InverseMat(T);
  h2Tinv:=h2*Tinv;
  if not h2Tinv in Sh1rec.horbit then 
    return(false);
  fi; 
  pos:=SinglePosition(Sh1rec.horbit, h2Tinv);
  tg:=Sh1rec.transporters[pos]*T;
  if TMTTmult(tg, GramS2)<> Sh1rec.GramS then beep(578517865); fi;
  if Sh1rec.h*tg<>h2 then beep(8121123); fi;
  return(tg);
end;


FoundShRecs:=[];

IsNewSh:=function(tK3rec)
  #
  local beep, isnewflag, nowSh, oldrec, tg, tginv, newspcons;
  #
  beep:=function(beepnumb)
    localbeep("IsNewSh", beepnumb); Error();
  end;
  #
  isnewflag:=true;
  nowSh:=[tK3rec.Gram, tK3rec.h];
  for oldrec in FoundShRecs do 
    tg:=IsIsomShs(oldrec, nowSh);
    if tg<>false then 
      isnewflag:=false;
      #
      tginv:=InverseMat(tg);
      newspcons:=(tK3rec.givenspcons)*tginv;
      if TMTTmult(newspcons, oldrec.Gram)<>tK3rec.tadj then beep(5817718); fi;
      Add(oldrec.spconss, newspcons);
      Add(oldrec.wgs, tK3rec.wg);
      #
      return();
    fi;
  od;
  #
  if isnewflag then 
    Add(FoundShRecs, GetIniData(tK3rec));
  fi;
  #
  return();
  #
end;

######