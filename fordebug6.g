

Gram:=[[2,2,2],[2,-2,0],[2,0,-2]];

h:=[1,0,0];

rats:=FindRatsOnK3(Gram, h, 2);

r1:=[0,1,0];
s1:=2*h-r1;
r2:=[0,0,1];
s2:=2*h-r2;
r3:=[1,1,-1];
s3:=2*h-r3;


A:=[[2,2,2,2],[2,-2,0,0], [2,0,-2,0], [2,0,0,-2]];
h:=MakeVectei(4, 1);

Gram:=MatMatToMat(
  [
   [ [[2]], [[2,2,2,2,2,2,2,2]] ], 
   [ TransposedMat([[2,2,2,2,2,2,2,2]]), (-2)*IdentityMat(8)]
   ]);

h:=MakeVectei(9, 1);

readdata("6TanConRecs");

Gram:=WeightedCompleteGraphToGram(3, [0,0,2]);
h:=MakeVectei(4, 1);
##########