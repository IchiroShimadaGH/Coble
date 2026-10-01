#Read("Task20260902Ub.g");


Read("Make_nRvRsRec.g");

thenRvRsrec:=Make_nRvRsRec(GramR, embR, GramL);

Printn("We prepare thenRvRsrec and discSrec");
Printn("thenRvRsrec");
for pos in [1..Length(thenRvRsrec.nRs)] do
  Printn(thenRvRsrec.nRs[pos], Length(thenRvRsrec.vRss[pos]));
od;

discSrec:=DiscriminantForm(GramS);

savedata(thenRvRsrec);
savedata(discSrec);
