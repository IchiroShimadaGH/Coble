#Read("BorcherdsCoble.g");


BorcherdsCoble:=function(datarec)
  local resultrec, beep;
  #
  beep:=function(beepnumb)
    localbeep("BorcherdsCoble", beepnumb); Error();
  end;
  #
  #
  gens:=[];
  chamrecs:=[];
  #
  #############
  # the main part
  #############
  for tchamrec in chamrecs do
  od;
  #
  resultrec:=rec(
    gens:=gens, chamrecs:=chamrecs
  );
  #
  return(resultrec);
end;

#####