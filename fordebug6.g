GetSpRats:=function(Gram, h)
  #
  local tinvol, ratss, sprats, spratsrec,  dones, tr, ttr, dd,
  rats;
  #
  tinvol:=RatInvolRecP2(Gram, h, h).invol;
  ratss:= NpriFindRatssOnK3(Gram, h, 2);
  dd:=0;
  spratsrec:=rec();
  for rats in ratss do
    dd:=dd+1;
    sprats:=[];
    dones:=[];
    for tr in rats  do
      if tr in dones then continue; fi;
      ttr:=tr*tinvol;
      if ttr<>tr then 
        Add(sprats, [tr, ttr]);
        Append(dones, [tr, ttr]);
      fi; 
    od;
    if dd=1 then spratsrec.slins:=sprats;
    elif dd=2 then spratsrec.scons:=sprats;
    else beep(78575);
    fi;
  od;
  spratsrec.nopss:=List([spratsrec.slins, spratsrec.scons], Length);
  return(spratsrec);
end;
