

#Read("Task20260920Ua.g");


Read("SplitConsTools.g");

tE8:=ADEHGram(["E8"]);
w:=[];

for ii in [1..8] do 
  for jj in [ii+1..8] do 
    if tE8[ii][jj]=0 then 
      Add(w, 2);
    elif tE8[ii][jj]=-1 then 
      Add(w, 0);
    else beep(69698169); 
    fi;
  od;
od;

trec:=WGraphToGramh(w, 8);
tGram:=trec.Gram;
th:=trec.h;

if IsGeom(tGram, th)<>true then buzz(5861); fi;
if NoOverlattice(tGram, th)<>true then buzz(2861); fi;
scons:=GetTotalSplcons(tGram, th);
Printn(Length(scons));


############