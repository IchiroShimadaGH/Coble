#Read("Task20260920Va.g");

tGram:=DiagonalMat(Concatenation([2], List([1..8], ii->-2)));
aa:=Concatenation([3], List([1..8], ii->-1));

addvv:=[0,1,1,1,1,0,0,0,0];
tGramdual:=InverseMat(tGram);
if addvv*tGramdual*addvv <>-2 then buzz(626351); fi;
H:=HermiteNormalFormIntegerMat(CopyAdd(tGram, addvv));
newbasisdual:=Filtered(H, tv->not IsZeroVect(tv));
newGram:=TMTTmult(newbasisdual, tGramdual);
newaa:=aa*tGram*InverseMat(newbasisdual);

rs:=AffESstd(newGram, newaa, 0, -2, true);
Printn(Length(rs));

for tr in rs do 
  Printn(tr*newbasisdual*tGramdual);
od;


addvv:=[1,1,1,1,1,1,0,0,0];
tGramdual:=InverseMat(tGram);
if addvv*tGramdual*addvv <>-2 then buzz(626351); fi;
H:=HermiteNormalFormIntegerMat(CopyAdd(tGram, addvv));
newbasisdual:=Filtered(H, tv->not IsZeroVect(tv));
newGram:=TMTTmult(newbasisdual, tGramdual);
newaa:=aa*tGram*InverseMat(newbasisdual);

rs:=AffESstd(newGram, newaa, 0, -2, true);
Printn(Length(rs));

for tr in rs do 
  Printn(tr*newbasisdual*tGramdual);
od;

addvv:=[0,1,1,1,1,1,1,1,1];
tGramdual:=InverseMat(tGram);
if addvv*tGramdual*addvv <>-4 then buzz(626351); fi;
H:=HermiteNormalFormIntegerMat(CopyAdd(tGram, addvv));
newbasisdual:=Filtered(H, tv->not IsZeroVect(tv));
newGram:=TMTTmult(newbasisdual, tGramdual);
newaa:=aa*tGram*InverseMat(newbasisdual);

rs:=AffESstd(newGram, newaa, 0, -2, true);
Printn(Length(rs));

for tr in rs do 
  Printn(tr*newbasisdual*tGramdual);
od;

addvv:=[1,1,0,0,0,0,0,0,0];
tGramdual:=InverseMat(tGram);
if addvv*tGramdual*addvv <>0 then buzz(626351); fi;
H:=HermiteNormalFormIntegerMat(CopyAdd(tGram, addvv));
newbasisdual:=Filtered(H, tv->not IsZeroVect(tv));
newGram:=TMTTmult(newbasisdual, tGramdual);
newaa:=aa*tGram*InverseMat(newbasisdual);

rs:=AffESstd(newGram, newaa, 0, -2, true);
Printn(Length(rs));

for tr in rs do 
  Printn(tr*newbasisdual*tGramdual);
od;

######