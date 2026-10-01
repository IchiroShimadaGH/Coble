
#Read("Task20260902Uc.g");

readdata("weyl0");
readdata("discSrec");
readdata("thenRvRsrec");

Borrec:=rec(
  GramS:=GramS, 
  GramR:=GramR, 
  GramL:=GramL, 
  embS:=embS, 
  embR:=embR, 
  projR:=projR, 
  projS:=projS,
  Rtype:=Rtype, 
  discS:=discS,
  discR:=discR,
  anisomdiscSR:=anisom ,
  GramSdual:=InverseMat(GramS), 
  GramRdual:=InverseMat(GramR),
  hS:=hS,
  hSL:=hS*embS, 
  aS:=aS, #an ample that is interior in cham 
  aSL:=aS*embS,
  aR:=aR, 
  aRL:=aR*embS,
  aL:=aL,
  weyl0:=weyl0, 
  discSrec:=discSrec
);

Borrec.embSdual:=Borrec.GramSdual*Borrec.embS;
Borrec.embRdual:=Borrec.GramRdual*Borrec.embR;
Borrec.projSdual:=Borrec.projS*Borrec.GramS;
Borrec.projRdual:=Borrec.projR*Borrec.GramR;

Borrec.nRvRsrec:=thenRvRsrec;

Printn("Borrec");
Printn("recnames", RecNames(Borrec));

savedata(Borrec);


