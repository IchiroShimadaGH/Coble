#
#Read("Task20260908Ua.g");



Read("NewOGLat.g");
Read("MassFormula.txt");

Rtype1:=List([1..9], ii->"A1");
Add(Rtype1, "D6");

GramR1:=ADEHGram(Rtype1);

themass1:=Mass(GramR1);

Rtype2:=List([1..7], ii->"A1");
Append(Rtype2, ["D4", "D4"]);

GramR2:=ADEHGram(Rtype2);

themass2:=Mass(GramR2);

if themass1<>themass2 then beep(77711); fi;

uvabR1:=uvabvect(DiscriminantForm(GramR1));

ogrec1:=RepeatNewOGLat(GramR1, 20, 1);;
ogrec2:=RepeatNewOGLat(GramR2, 20, 1);;

Printn(ogrec1.size);
Printn(ogrec2.size);

ogrecs:=[ogrec1, ogrec2];

tmass:=(1/ogrec1.size)+(1/ogrec2.size);


specialRepeatNewOGLat:=function(GramL, basisrectrial, giveupsec)
  local counter, basisrec, OGrec, tab, nn, tU;
  counter:=0;
  tab:=10;
  nn:=Length(GramL);
  while true do
    counter:=counter+1;
    tU:=RandomUnimodMat(nn);
    basisrec:=BasisRec(TMTTmult(tU, GramL),  basisrectrial+tab*counter);
    OGrec:=NewOGLat(basisrec, giveupsec+counter);
    if OGrec<>fail then 
      OGrec.basisrectrial:= basisrectrial+counter;
      OGrec.giveupsec:=giveupsec+counter;
      return(OGrec);
    else 
      Printn("____repeat", counter);
    fi;
  od;
end;

tempsavename:="togrecs2";


while tmass<themass1 do 
  rogrec:=Random(ogrecs);
  tGramL:=RandomKneserpNbrL(rogrec.Gram, 3);
  if not isisomuvab(uvabR1, uvabvect(DiscriminantForm(tGramL))) then 
    beep(76171781);
  fi;
  if Mass(tGramL)<>themass1 then 
    beep(44471781);
  fi;
  tbasiscrec:=BasisRec(tGramL, 20);
  isnewflag:=true;
  for oldogrec in ogrecs do 
    if IsIsomBasisRecs(oldogrec.basisrec, tbasiscrec)<>false then 
      isnewflag:=false;
      Printn("____not new");
      break; #from for oldogrec in ogrecs do 
    fi;
  od;
  if isnewflag then 
    newogrec:=specialRepeatNewOGLat(tGramL, 10, 100);;  #10 times LLLtrial: initial giveup 100 sec
    Add(ogrecs, newogrec);
    tmass:=tmass+1/(newogrec.size);
    Printn("new", newogrec.size);
    savedataas(ogrecs, tempsavename);
  fi;
od;
#