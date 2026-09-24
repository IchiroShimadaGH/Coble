#Read("ConcatrandomPro.g");

names := DirectoryContents("gapsaves");;
files := Filtered(names, name -> PositionSublist(name, "saverandomPro") = 1);


total:=[];
for tf in files do 
  ttf := tf{[5 .. Length(tf)-4]};
  readdata(ttf);
  Append(total, EvalString(ttf));
od;


Printn(Length(total));

nopss:=List(total, xx->xx[1]);
Sort(nopss);
nopss:=Reversed(nopss);


Printn(List([1..20], ii->nopss[ii]));
