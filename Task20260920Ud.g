

#Read("Task20260920Ud.g");


Read("SplitConsTools.g");

adjE8:=[ [ -2, 2, 2, 0, 2, 2, 2, 2 ], 
[ 2, -2, 0, 2, 2, 2, 2, 2 ], 
[ 2, 0, -2, 0, 2, 2, 2, 2 ], 
[ 0, 2, 0, -2, 0, 2, 2, 2 ], 
[ 2, 2, 2, 0, -2, 0, 2, 2 ], 
[ 2, 2, 2, 2, 0, -2, 0, 2 ], 
[ 2, 2, 2, 2, 2, 0, -2, 0 ], 
[ 2, 2, 2, 2, 2, 2, 0, -2 ] ];

specialnewscons3:=[];

addnumb:=3;
counter:=0;
while true do 
  counter:=counter+1;
  if counter mod 1000000=0 then Printn("counter", counter); fi;
  tA:=RandomMatFromL(addnumb, 8, [0,1,2,3,4]);
  ttA:=TransposedMat(tA);
  tB:=RandomMatFromL(addnumb, addnumb, [0,1,2,3,4]);
  for ii in [1..addnumb] do tB[ii][ii]:=-2; od;
  newadj:=MatMatToMat([[adjE8, ttA], [tA, tB]]);
  trec:=WGraphToGramh(newadj, 8+addnumb);
  tGram:=trec.Gram;
  th:=trec.h;
  isgeom:=IsGeom(tGram, th);
  if isgeom<>true then  
    # if isgeom<>[ false, 2 ]  then 
    #   Printn(counter, "not geom", isgeom);
    # fi;
    continue;
  fi;
  if NoOverlattice(tGram, th)<>true then  
    Printn("there exists an overlattice");
  fi;
  scons:=GetTotalSplcons(tGram, th);
  Printn(Length(scons), tA, tB);
  if Length(scons)>120 then 
    Add(specialnewscons3, [tA, tB]);
    savedata(specialnewscons3); 
    Printn("___", counter, Length(specialnewscons3));
  fi;
  if Length(scons)<120 then beep(761476); fi;
od;


############