#Read("Task20260901Ua.g");



aa:=List([1..10], ii->-2);
Add(aa, 2, 1);
Printn(aa);

GramS:=DiagonalMat(aa);
discS:=DiscriminantForm(GramS);
discmS:=DiscriminantForm(-GramS);

uvabvectS:=uvabvect(discS);
uvabvectmS:=uvabvect(discmS);

#Printn(isisomuvab(uvabvectS, uvabvectmS));




ADs:=["A1"];
ADranks:=[1];
for hrr in [2..7] do 
  rr:=2*hrr;
  Drr:=Concatenation("D", String(rr));
  Add(ADs, Drr);
  Add(ADranks, rr);
od;


candidates:=[];

task:=function(ttype, trank)
  local tGram, newrank, muvab, mdisc;
  if trank=15 then 
    if ADEHrank(ttype)<>15 then beep(116986981); fi;
    tGram:=-ADEHGram(ttype);#newgative-definite
    mdisc:=DiscriminantForm(-tGram); # MINUS q
    muvab:=uvabvect(mdisc);
    if isisomuvab(uvabvectS, muvab) then 
      Printn("found", ttype);
      sttype:=Collected(ttype);
      AddSet(candidates, sttype);
    fi;
  else 
    for tpos in [1..Length(ADs)] do
      Add(ttype, ADs[tpos]);
      newrank:=ADranks[tpos];
      if trank+newrank<=15 then 
        task(ttype, trank+newrank);
      fi;
      Remove(ttype);
    od;
  fi;
end;

task([], 0);

Printn(candidates);

#####