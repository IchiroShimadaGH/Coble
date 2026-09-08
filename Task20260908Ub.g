#
#Read("Task20260908Ub.g");



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

IsIsomBasisRecsADE:=function(basisrecA, basisrecB)
  if basisrecA.ADEtype<>basisrecB.ADEtype then return(false); fi;
  return(IsIsomBasisRecs(basisrecA, basisrecB));
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


tempsavename:="basisrecs2";

basisrecs:=[basisrecR1, basisrecR2];

while true do 
  rbasisrec:=Random(basisrecs);
  tGramL:=RandomKneserpNbrL(rbasisrec.Gram, 3);
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
    if IsIsomBasisRecsADE(oldbasisrec, tbasisrec)<>false then 
      isnewflag:=false;
      Printn("____not new");
      break; #from for oldogrec in ogrecs do 
    fi;
  od;
  if isnewflag then 
    Add(basisrecs, tbasisrec);
    ttype:=tbasisrec.ADEtype;
    Printn("new", Length(basisrecs), Collected(ttype), ADEHrank(ttype));
    savedataas(basisrecs, tempsavename);
  fi;
od;
#