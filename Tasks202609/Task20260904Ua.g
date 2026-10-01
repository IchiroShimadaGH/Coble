#Read("Task20260904Ua.g");

Read("Adj.g");
readdata("Borrec");
readdata("nSaSvSss0");

walls0:=Concatenation(List(nSaSvSss0, xx->List(xx.vSliftss, yy->yy[1])));

if Set(walls0, tv->IsWalldefv(Borrec, tv))<>[true] then beep(818181); fi;

weyl0:=Borrec.weyl0;

adjrecs0:=[];

for tv in walls0 do 
  adjrec:=AdjWeyl(Borrec, weyl0, tv);
  Add(adjrecs0, adjrec);
od;

#CheckAdjrec:=function(Borrec, adjrec, trialnumb)

pos:=0;
for tadjrec in adjrecs0 do 
  pos:=pos+1;
  CheckAdjrec(Borrec, tadjrec, 5);
  #Printn(pos, Length(adjrecs0));
od;

Printn("adjrecs0 is checked");
savedatap(adjrecs0);

####################