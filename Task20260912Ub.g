#Read("Task20260912Ub.g");

hh2:=[3,-1,-1,-1,-1,-1,-1,-1, -1];

if hh2*GramS*hh2<>2 then buzz(112134); fi;

readdata("240rats");

hh2dual:=hh2*GramS;

if Set(240rats*hh2dual)<>[2] then buzz(771621); fi;


#function( SQ, vec, reldeg, norm, isES ) ... end
if AffESstd(GramS, hh2, 0, -2, true)<>[] then buzz(51232); fi;

aa:=hh2;
aadual:=aa*GramS;

if SeparatingVects(GramS, hh, aa, -2)<>[] then buzz(11422); fi;



invol2:=RatInvolRecP2(GramS, aa, aa).invol;;


hhdual:=hh*GramS;
Printn(Collected(240rats*hhdual));
Set(240rats, tr->tr*GramS*(tr*invol2));


dones:=[];
orbs:=[];

for tr in 240rats do 
  if not tr in dones then 
    orb:=[tr, tr*invol2];
    Add(orbs, orb);
    Append(dones, orb);
  fi;
od;

6tancons:=List(orbs);

intpat6tcs:=function(6tc1, 6tc2)
  local 6tc1dual, pm, tv1, tv2, type;
  6tc1dual:= 6tc1*GramS;
  pm:=List(6tc1dual, tv1->List(6tc2, tv2->tv1*tv2));
  if pm=[[4,0],[0,4]] then  type:=0; 
  elif pm=[[3,1],[1,3]] then type:=1; 
  elif pm=[[2,2],[2,2]] then type:=2;
  elif pm=[[1,3],[3,1]] then type:=1;
  elif pm=[[0,4],[4,0]] then type:=0;
  else 
    if not IsEqualSet(6tc1, 6tc2) then buzz(61615); 
    else 
      if pm<>[[-2, 6], [6, -2]] and  pm<>[[6,-2], [-2, 6]] then buzz(33315); fi;
      type:=-2;
    fi;
  fi;
  return(type);
end;


bigPM:=List(6tancons, tc1->List(6tancons, tc2->intpat6tcs(tc1, tc2)));


if bigPM<>TransposedMat(bigPM) then beep(991919); fi;
Set( bigPM, Collected)

#####