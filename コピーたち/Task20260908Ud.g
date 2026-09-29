#
#Read("Task20260908Ud.g");



Read("NewOGLat.g");
Read("MassFormula.txt");

AddADEtype:=function(basisrec)
  if 2 in basisrec.nrmsset then 
    basisrec.ADEtype:=GetRootType(basisrec.Gram, basisrec.vss[1]);
    Sort(basisrec.ADEtype);
  else 
    basisrec.ADEtype:=[];
  fi;
  return();
end;




Rtype1:=List([1..9], ii->"A1");
Add(Rtype1, "D6");
GramR1:=ADEHGram(Rtype1);
basisrecR1:=BasisRec(GramR1, 100);
AddADEtype(basisrecR1);

Rtype2:=List([1..7], ii->"A1");
Append(Rtype2, ["D4", "D4"]);
GramR2:=ADEHGram(Rtype2);
basisrecR2:=BasisRec(GramR2, 100);
AddADEtype(basisrecR2);


themass:=Mass(GramR1);
uvabR:=uvabvect(DiscriminantForm(GramR1));


EasyIsIsomBasisRecs:=function(basisrecA, basisrecB)
  #
  local GramL1, GramL2, nn, vss2, vss2dual, getvs, 
  getvsdual, basis1, fingerprints1, basis1inv, totalflag, newGram1, 
  thetg, extend, cleng, ii, revflag, basisrec1, basisrec2;
  #
  if basisrecA.ADEtype<>basisrecB.ADEtype then return(false); fi;
  if basisrecA.minflag and basisrecB.minflag then 
    if basisrecA.nrmsset<>basisrecB.nrmsset then return(false); fi;
    if basisrecA.nopss<>basisrecB.nopss then return(false); fi;
    revflag:=false;
    basisrec1:=basisrecA;
    basisrec2:=basisrecB;
  else 
    if IsSubset(basisrecA.nrmsset, basisrecB.nrmsset) then 
      cleng:=Length(basisrecB.nrmsset);
      revflag:=true;
      basisrec1:=basisrecB;
      basisrec2:=basisrecA;
    elif IsSubset(basisrecB.nrmsset, basisrecA.nrmsset) then 
      cleng:=Length(basisrecA.nrmsset);
      revflag:=false;
      basisrec1:=basisrecA;
      basisrec2:=basisrecB;
    else 
      return(false);
    fi;
    for ii in [1..cleng] do 
      if basisrecA.nopss[ii]<>basisrecB.nopss[ii] then 
        return(false);
      fi;
    od;
  fi;
  #
  #
  GramL1:=basisrec1.Gram;
  GramL2:=basisrec2.Gram;
  nn:=Length(GramL1);
  if nn<>Length(GramL2) then return(false); fi;
  if DeterminantIntMat(GramL1)<>DeterminantIntMat(GramL2) then return(false); fi;
  #
  return(true);
  #
end;

tempsavename:="easybasisrecs5b";

basisrecs:=[basisrecR1, basisrecR2];

notnewcounter:=0;

while true do 
  rbasisrec:=Random(basisrecs);
  tGramL:=RandomKneserpNbrL(rbasisrec.Gram, 5);
  if not isisomuvab(uvabR, uvabvect(DiscriminantForm(tGramL))) then 
    beep(76171781);
  fi;
  if Mass(tGramL)<>themass then 
    beep(44471781);
  fi;
  tbasisrec:=BasisRec(tGramL, 100);
  AddADEtype(tbasisrec);
  isnewflag:=true;
  for oldbasisrec in basisrecs do 
    if EasyIsIsomBasisRecs(oldbasisrec, tbasisrec)<>false then 
      isnewflag:=false;
      #Printn("____not new");
      notnewcounter:=notnewcounter+1;
      break; #from for oldogrec in ogrecs do 
    fi;
  od;
  if isnewflag then 
    notnewcounter:=0;
    Add(basisrecs, tbasisrec);
    ttype:=tbasisrec.ADEtype;
    Printn("new", Length(basisrecs), Collected(ttype), ADEHrank(ttype));
    savedataas(basisrecs, tempsavename);
  fi;
  if notnewcounter>0 and notnewcounter mod 100=0 then Printn("____not new", notnewcounter); fi;
od;
#