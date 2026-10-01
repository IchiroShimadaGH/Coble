#Read("Task20260901Uc.g");



mindeg:=infinity;
LampleSs:=[];

while Length(LampleSs)<10000 do
  rv:=RandomVectFromL(11, [-3, -2, -1,-1,0,1,  1, 2, 3]);
  kk:=0;
  cvRss:=[];
  while true do 
    kk:=kk+1;
    bS:=kk*aS+rv;
    if bS*GramS*bS>0 then
      bSL:=bS*embS;
      tvs:=AffESstd(GramL,  bSL, 0, -2, true);;
      cvRs:=Collected(List(tvs*projR, tx->tx*GramR*tx));
      newdeg:=bS*GramS*bS;
      Printn(newdeg,  cvRs, Length(LampleSs), mindeg);
      if cvRs=[[-2, nrR]] then 
        if newdeg<mindeg then 
          LampleS:=bS;
          mindeg:=newdeg;
          savedata(LampleS);
        fi;
        mindeg:=Minimum(newdeg, mindeg);
        Add(LampleSs, bS); 
        savedata(LampleSs);
        break; 
      else 
        if cvRs in cvRss then break; 
        else 
          Add(cvRss, cvRs);
        fi;
      fi;
    fi;
  od;
od;

Printn(List(LampleSs, tv->tv*GramS*tv));




