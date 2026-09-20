

#Read("Task20260920Vc.g");


Read("SplitConsTools.g");



specialnewscons11:=[];

kk:=11;
counter:=0;
while true do 
  counter:=counter+1;
  if counter mod 10000=0 then Printn("counter", counter); fi;
  newadj:=RandomSymMatFromL(11, [0,1,2,3,4]);
  for ii in [1..kk] do newadj[ii][ii]:=-2; od;
  trec:=WGraphToGramh(newadj, kk);
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
  Printn(Length(scons));
  if Length(scons)>120 then 
    Add(specialnewscons11, newadj);
    savedata(specialnewscons11); 
    Printn("___", counter, Length(specialnewscons11));
  fi;
  
od;


############