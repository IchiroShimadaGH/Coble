

#Read("Task20260920Uc.g");


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

counter:=0;
while true do 
  counter:=counter+1;
  if counter mod 10000=0 then Printn("counter", counter); fi;
  tv1:=RandomVectFromL(8, [0,1,2,3,4]);
  tv2:=RandomVectFromL(8, [0,1,2,3,4]);
  aa:=Random([0,1,2,3,4]);
  tA:=[tv1, tv2];
  ttA:=TransposedMat(tA);
  tB:= [[-2, aa], [aa, -2]];
  newadj:=MatMatToMat([[adjE8, ttA], [tA, tB]]);
  trec:=WGraphToGramh(newadj, 10);
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
  Printn(Length(scons), tA);
  if Length(scons)>120 then 
    Add(specialnewscons2, [tA, tB]);
    savedata(specialnewscons2); 
    Printn("___", counter, Length(specialnewscons2));
  fi;
  if Length(scons)<120 then beep(761476); fi;
od;


############