#Read("Task20260904Ue.g");

Read("GeneralStabChain.g");
Read("Make_nSaSvSss.g");
Read("AutCham.g");
Read("Walls.g");
Read("Adj.g");

readdata("Borrec");

GramS:=Borrec.GramS;
projS:=Borrec.projS;
GramSdual:=InverseMat(GramS);

ChamNumDatas2:=[];

tweyls:=[Borrec.weyl0];
tweylSs:=[Borrec.weyl0*projS];

pos:=0;

for tweyl in tweyls do 
  pos:=pos+1;
  tnSaSvSss:=Make_nSaSvSss(Borrec, tweyl);
  #check_nSaSvSss(Borrec,  tweyl, tnSaSvSss);
  AddSet(ChamNumDatas2, ChamNumData(tnSaSvSss));
  twalls:=Concatenation(List(tnSaSvSss, xx->List(xx.vSliftss, yy->yy[1])));
  twallslf:=Set(twalls*GramS, GetPrimitiveIntVect);
  wallslf:=BorcherdsSelectWalls(twallslf);
  if Length(twallslf)<>Length(twallslf) then Printn("LP does a job!"); fi;
  walls:=wallslf*GramSdual;
  for twall in walls do  
    adjrec:=AdjWeyl(Borrec, tweyl, twall);
    #CheckAdjrec(Borrec, adjrec, 3);
    newweyl:=adjrec.newweyl;
    newweylS:=newweyl*projS;
    if not newweylS in tweylSs then 
      AddSet(tweylSs, newweylS);
      Add(tweyls, newweyl);
    fi;
  od;
  Printn(pos, Length(walls), Length(tweyls), Length(ChamNumDatas2));
  savedata(ChamNumDatas2);
od;


####################