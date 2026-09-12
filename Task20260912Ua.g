#Read("Task20260912Ua.g");


buzz:=function(numb)
  Printn("Error buzz", numb); Error();

end;

dd:=[2,-2,-2,-2,-2,-2,-2,-2,-2];;
if Collected(dd)<>[[-2,8], [2,1]] then buzz(3581); fi;

GramS:=DiagonalMat(dd);

hh:=MakeVectei(9, 1);
exceps:=List([2..9], ii->MakeVectei(9, ii));


FindRatsOnK3Excep:=function(GramSX, pol, exceps, uptodeg)
  local rats, extradeg, ratsdual, task, tdeg;
  #
  extradeg:=0;
  rats:=List(exceps);
  ratsdual:=rats*GramSX;
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
  tdeg:=0;
  #
  while tdeg<uptodeg do
    tdeg:=tdeg+1;
    AffES(GramSX, pol, tdeg, -2, true, task);
    Printn("tdeg", tdeg, Length(rats));
    if Length(rats)=240 then savedataas(rats, "240rats");fi;
    if Length(rats)>240 then buzz(613211); fi;
  od;
  #
  return(rats);
end;


FindRatsOnK3Excep(GramS, hh, exceps, 20);




######