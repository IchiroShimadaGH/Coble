
#Read("Task20260927Ua.g");


Read("MinimalWgs.g");

EkkRecs:=[[]];

counter:=0;
while true do
  for kk in [2..10] do 
    EkkRec:=GetEkkRec(kk);
    rwg:=RandomVectFromL(EkkRec.leng, [-3..3]);
    tau1:=Random(SymmetricGroup(kk));
    tau2:=Random(SymmetricGroup(kk));
    #
    tau12:=tau1*tau2;
    tau21:=tau2*tau1;
    tau1inv:=Inverse(tau1);
    tau2inv:=Inverse(tau2);
    tau1inv2:=tau1inv*tau2;
    tau12inv:=tau1*tau2inv;
    tau1inv2inv:=tau1inv*tau2inv;
    #
    ltau1:=TauToLargetau(kk, EkkRec.Ekk, tau1);
    ltau2:=TauToLargetau(kk, EkkRec.Ekk, tau2);
    ltau12:=TauToLargetau(kk, EkkRec.Ekk, tau12);
    ltau21:=TauToLargetau(kk, EkkRec.Ekk, tau21);
    ltau1inv:=TauToLargetau(kk, EkkRec.Ekk, tau1inv);
    ltau2inv:=TauToLargetau(kk, EkkRec.Ekk, tau2inv);
    ltau1inv2:=TauToLargetau(kk, EkkRec.Ekk,  tau1inv2);
    ltau12inv:=TauToLargetau(kk, EkkRec.Ekk,  tau12inv);
    ltau1inv2inv:=TauToLargetau(kk, EkkRec.Ekk,  tau1inv2inv);
    #
    if ltau12<>ltau1*ltau2 then buzz(7666476); fi;
    if ltau21<>ltau2*ltau1 then buzz(3666476); fi;
    if ltau1inv<>Inverse(ltau1) then buzz(761476); fi;
    if ltau2inv<>Inverse(ltau2) then buzz(733476); fi;
    if ltau12inv<>ltau1*ltau2inv then buzz(3366476); fi;
    if ltau1inv2<>ltau1inv*ltau2 then buzz(5566476); fi;
    if ltau1inv2inv<>ltau1inv*ltau2inv then buzz(996476); fi;
    #
    if ActionTauWg(ltau12, rwg)<>ActionTauWg(ltau1, ActionTauWg(ltau2, rwg)) then buzz(9964336); fi;
    if ActionTauWg(ltau21, rwg)<>ActionTauWg(ltau2, ActionTauWg(ltau1, rwg)) then buzz(9964336); fi;
    #
    if rwg<>ActionTauWg(ltau1, ActionTauWg(ltau1inv, rwg)) then buzz(9964336); fi;
    if rwg<>ActionTauWg(ltau2, ActionTauWg(ltau2inv, rwg)) then buzz(997736); fi;
    if rwg<>ActionTauWg(ltau1inv, ActionTauWg(ltau1, rwg)) then buzz(922336); fi;
    if rwg<>ActionTauWg(ltau2inv, ActionTauWg(ltau2, rwg)) then buzz(988336); fi;
    #
    if ActionTauWg(ltau12inv, rwg)<>ActionTauWg(ltau1, ActionTauWg(ltau2inv, rwg)) then
       buzz(99634336); 
    fi;
    if ActionTauWg(ltau1inv2, rwg)<>ActionTauWg(ltau1inv, ActionTauWg(ltau2, rwg)) then
       buzz(9963334336); 
    fi;
    if ActionTauWg(ltau1inv2inv, rwg)<>ActionTauWg(ltau1inv, ActionTauWg(ltau2inv, rwg)) then
       buzz(995563334336); 
    fi;
  od;
  counter:=counter+1;
  if counter mod 1000=0 then Printn(counter); fi;
od;







