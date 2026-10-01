#Read("Task20260917Ub.g");

Read("Wk.g");



for tt in [1..1000000] do 
  for kk in [2..5] do 
    if kk=2 then Gorbs:=GorbsW2; 
    elif kk=3 then Gorbs:=GorbsW3;
    elif kk=4 then Gorbs:=GorbsW4;
    elif kk=5 then Gorbs:=GorbsW5;
    fi;
    #
    cc:=2^(kk-2);
    leng:=kk*(kk-1)/2;
    for ss in [1..cc] do 
      tw:=RandomVectFromL(leng, [0,1,2,3,4]);
      mtw:=MinRepWk(kk, tw);
      pos:=SinglePosition(Gorbs, mtw);
      if OrbitWk(kk, tw)<>OrbitWk(kk, mtw) then buzz(7171); fi;
    od;
  od;
  if tt mod 100=0 then Printn(tt);fi;
od;

# GorbsW2:=MinRepsListWk(2);
# Printn("k=2", Length(GorbsW2));
# GorbsW3:=MinRepsListWk(3);
# Printn("k=3", Length(GorbsW3));
# GorbsW4:=MinRepsListWk(4);
# Printn("k=4", Length(GorbsW4));
# GorbsW5:=MinRepsListWk(5);
# Printn("k=5", Length(GorbsW5));

# savedata(GorbsW2);
# savedata(GorbsW3);
# savedata(GorbsW4);
# savedata(GorbsW5);

# gap> Read("Task20260917Ua.g");
# k=2 3 
# k=3 14 
# k=4 161 
# k=5 6595


#####