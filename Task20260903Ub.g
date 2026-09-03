#Read("Task20260903Ub.g");

# readdata("Borrec");
# readdata("nSaSvSss0");
# readdata("autchamrec0");

#Read("BorcherdsSelectWalls.g");
Read("Walls.g");

twalls0:=Concatenation(List(nSaSvSss0, xx->List(xx.vSliftss, yy->yy[1])));

if Length(twalls0)<>41 then beep(9919191); fi;

GramS:=Borrec.GramS;

nnss:=[];
nops:=0;
pos1:=0;
for tv1 in twalls0  do
  pos1:=pos1+1;
  tv1dual:=tv1*GramS;
  pos2:=0;
  for tv2 in twalls0  do
    pos2:=pos2+1;
    if pos1<pos2 then 
      aa:=tv1dual*tv2;
      if aa>0 then 
        nops:=nops+1; 
        nns:=[tv1dual*tv1, tv2*GramS*tv2];
        Add(nnss, nns);
      fi;
    fi;
  od;
od;

Printn(nops);
Printn(Collected(nnss));

Printn("__________________________");

twalls0lf:=Set(twalls0*GramS);;

if Set(twalls0lf, IsIntVect)<>[true] then beep(51515); fi;
if  Set(twalls0lf, IsPrimitiveVect)<>[true] then beep(22515); fi;

if twalls0lf<>Set(twalls0lf, lf->GetPrimitiveIntVect(lf)) then beep(9991111); fi;


walls0lf:=BorcherdsSelectWalls(twalls0lf);

if walls0lf<>twalls0lf then beep(919911); fi;

wallsresults:=Make_Walls(Borrec, autchamrec0.periodorbs);

if Set(wallsresults, xx->xx[1]) <>[true] then beep(887711); fi;

Printn("Every element in walls0 is really a wall.");
Printn("Calculated by BorcherdsSelectWalls and by Make_Walls.");

w0S:=Borrec.weyl0*Borrec.projS;

if Set(walls0lf*w0S, xx->xx>0)<>[true] then beep(51112); fi;
Printn("__________________________");
Printn("w0S:=Borrec.weyl0*Borrec.projS is an interior point of cham 0.");



#####