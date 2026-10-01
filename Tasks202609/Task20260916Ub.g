#Read("Task20260916Ub.g");

Read("WeightedCompleteGraphClasses.g");



counter:=0;
for kk in [2..4] do 
  tgras:=WeightedCompleteGraphClasses(kk);
  for tgra in tgras do 
    counter:=counter+1;
    tGram:=WeightedCompleteGraphToGram(kk, tgra);
    if Rank(tGram)=kk+1 then 
      #Printn("nondeg", tgra, counter); 
    else 
      Printn("deg", tgra, counter,  kk+1-Rank(tGram)); 
    fi;
  od;
od;


#####