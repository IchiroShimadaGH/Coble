

#Read("Task20260920Vb.g");


Read("SplitConsTools.g");

adjE8:=[ [ -2, 2, 2, 0, 2, 2, 2, 2 ], 
[ 2, -2, 0, 2, 2, 2, 2, 2 ], 
[ 2, 0, -2, 0, 2, 2, 2, 2 ], 
[ 0, 2, 0, -2, 0, 2, 2, 2 ], 
[ 2, 2, 2, 0, -2, 0, 2, 2 ], 
[ 2, 2, 2, 2, 0, -2, 0, 2 ], 
[ 2, 2, 2, 2, 2, 0, -2, 0 ], 
[ 2, 2, 2, 2, 2, 2, 0, -2 ] ];

specialnewscons2:=[];

addnumb:=2;
counter:=0;
while true do 
  counter:=counter+1;
  if counter mod 10000=0 then Printn("counter", counter); fi;
  tA:=RandomMatFromL(addnumb, 8, [0,1,2,3,4]);
  ttA:=TransposedMat(tA);
  tB:=RandomMatFromL(addnumb, addnumb, [0,1,2,3,4]);
  for ii in [1..addnumb] do 
    tB[ii][ii]:=-2; 
    for jj in [ii+1..addnumb] do 
      tB[jj][ii]:=tB[ii][jj];
    od;
  od;
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
  Printn(counter, "geom");
  scons:=GetTotalSplcons(tGram, th);
  Printn(Length(scons), tA, tB);
  if Length(scons)>120 then 
    Add(specialnewscons2, [tA, tB]);
    savedata(specialnewscons2); 
    Printn("___", counter, Length(specialnewscons2));
  fi;
  if Length(scons)<120 then beep(761476); fi;
od;


############