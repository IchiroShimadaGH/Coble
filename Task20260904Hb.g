#Read("Task20260904Hb.g");


Stype:=["H2"];
for ii in [1..10] do 
  Add(Stype, "A1");
od;

discS:=NegADEHdiscrec(Stype);

diagsS:=[2];
for ii in [1..10] do 
  Add(Stype, -2);
od;

GramS:=DiagonalMats(diagsS);

Rtype:=["D6"];
for ii in [1..9] do 
  Add(Rtype, "A1");
od;

discR:=NegADEHdiscrec(Rtype);
mdiscR:=PosADEHdiscrec(Rtype);

mGramR:=ADEHGram(Rtype);
GramR:=-mGramR;

# discS is a +10b 
# discR is d6+9a

tA:=[ [ 1, 0 ], [ 1, 1 ] ]; #from 2b to posd6
tB:=[[0,0,1], [1,0,0], [1,1,0]];# from a+2b  to posd6+a 
tC:=[[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]; #from 4b to 4a 

isomSR:=DiagonalMats([tB, tC, tC]);
isomRS:=InverseMat(isomSR) mod 2;

if CopyNormalDiscf(TMTTmult(isomSR, mdiscR.discf))<>discS.discf then 
  beep(919191);
fi;

if CopyNormalDiscf(TMTTmult(isomRS, discS.discf))<>mdiscR.discf then 
  beep(922291);
fi;



OGRrec:=OGLat(mGramR);;

if OGRrec.order<>WeylGroupOrder(Rtype)*Factorial(9)*2 then 
  beep(1989191);
fi;

# gap> OGRrec.order;
# 8561413324800
# gap> uvabOqLOrder([0,0,1,10]);
# 46998591897600
# gap> 46998591897600/Factorial(9)/2;
# 64757760
uvabOqLOrder([0,0,0,2]);




#####