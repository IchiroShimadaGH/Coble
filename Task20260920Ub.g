

#Read("Task20260920Ub.g");


Read("SplitConsTools.g");

adjE8:=[ [ -2, 2, 2, 0, 2, 2, 2, 2 ], 
[ 2, -2, 0, 2, 2, 2, 2, 2 ], 
[ 2, 0, -2, 0, 2, 2, 2, 2 ], 
[ 0, 2, 0, -2, 0, 2, 2, 2 ], 
[ 2, 2, 2, 0, -2, 0, 2, 2 ], 
[ 2, 2, 2, 2, 0, -2, 0, 2 ], 
[ 2, 2, 2, 2, 2, 0, -2, 0 ], 
[ 2, 2, 2, 2, 2, 2, 0, -2 ] ];

specialnewscons:=[];

counter:=0;
while true do 
  counter:=counter+1;
  tv:=[RandomVectFromL(8, [0,1,2,3,4])];
  ttv:=TransposedMat(tv);
  newadj:=MatMatToMat([[adjE8, ttv], [tv, [[-2]] ]]);
  trec:=WGraphToGramh(newadj, 9);
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
  Printn(Length(scons), tv);
  if Length(scons)>120 then 
    Add(specialnewscons, tv);
  fi;
  if Length(scons)<120 then beep(761476); fi;
  if Length(specialnewscons) mod 100=0 then 
    Printn("___", counter, Length(specialnewscons));
    savedata(specialnewscons); 
  fi;
od;


############