#
#Read("Task20260908Va.g");



Read("NewOGLat.g");
Read("MassFormula.txt");




Rtype1:=List([1..9], ii->"A1");
Add(Rtype1, "D6");
GramR1:=ADEHGram(Rtype1);
basisrecR1:=BasisRec(GramR1, 100);

themass:=Mass(GramR1);
uvabR:=uvabvect(DiscriminantForm(GramR1));

readdata("easybasisrecs2");
if Length(easybasisrecs2)<>25 then beep(9191919);
fi;

Fs:=List(easybasisrecs2, xx->xx.Gram);
savedata(Fs);

sizes:=
[ 8561413324800, 1712282664960, 5096079360, 743178240, 28665446400, 1698693120,
298598400, 127401984, 452984832, 99532800, 2090188800, 743178240, 707788800,
3483648000, 67947724800, 16307453952, 321052999680, 9555148800, 133772083200,
31708938240, 535088332800, 67947724800, 570760888320, 2853804441600,
64210599936000 ];


#