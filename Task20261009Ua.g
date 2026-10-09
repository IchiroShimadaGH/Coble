
 #Read("Task20261009Ua.g");


#readdata("tempShrecs");

Length(tempShrecs);
#readdata("FoundShRecsOverLats");


#gap> Set(tempShrecs, RecNames);
#[ [ "h", "rank", "GramS", "GramP", "AutShgens", "AutShsize", 
#"Nikulin", "basisrecP", "detS", "horbit", "ogrecP", "sprats", "transporters" ] ]
#

Collected(List(tempShrecs, xx->xx.rank));
#[ [ 3, 3 ], [ 4, 9 ], [ 5, 45 ], [ 6, 326 ], [ 7, 5192 ] ] 

Collected(List(tempShrecs, xx->xx.sprats));

sconstypesTovect:=function(sconstypes)
  local zfttoddvect, xx;
  zfttoddvect:=[0,0,0];
  for xx in sconstypes do 
    if xx[1]="zf" then zfttoddvect[1]:=xx[2]; 
    elif xx[1]="tt" then zfttoddvect[2]:=xx[2];
    elif xx[1]="odd" then zfttoddvect[3]:=xx[2];
    else buzz(476476); 
    fi;
  od;
  return(zfttoddvect);
end;


tabledatas:=[];


for trec in tempShrecs do 
  trank:=trec.rank;
  th:=trec.h;
  if th*trec.GramS*th<>2 then beep(626872); fi;
  gens:=[th];
  for tpcs in trec.sprats.scons do 
    if Sum(tpcs)<>2*th then buzz(7162512); fi;
    Add(gens, tpcs[1]);
  od;
  if Rank(gens)<>trank then beep(144276); fi;
  coker:=CokerTorsion(gens);
  if coker=[] then extdeg:=1;
  else extdeg:=Product(coker);
  fi;
  typevect:=sconstypesTovect(trec.sprats.sconstypes);
  discS:=AbsInt(DeterminantIntMat(trec.GramS));
  tabledata:=[trec.rank, discS, extdeg, trec.sprats.nopss[2], trec.sprats.nopss[1], typevect];
  if extdeg>1 then Printn(tabledata); fi;
  Add(tabledatas, tabledata);
od;

Sort(tabledatas);

savedata(tabledatas);

#tabledatas2 was made by "temp2Shrecs"
#if tabledatas<>tabledatas2 then beep(55812); fi;



#######################