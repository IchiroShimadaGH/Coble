#Read("Task20260905Ub.g");

diags:=[1];
for ii in [1..10] do Add(diags, -1); od;

GramS11:=DiagonalMat(diags);

talphadual:=[];#ANNW
talpha:=[];#ANNW

sroots:=[];
refls:=[];

idmat11:=IdentityMat(11);
zv:=List([1..11], kk->0);

alpha0:=List(zv);
alpha0[1]:=1; alpha0[2]:=-1; alpha0[3]:=-1; alpha0[4]:=-1;
alpha0dual:=alpha0*GramS11;
if alpha0dual*alpha0<>-2 then beep(919191); fi;
R0:=List(idmat11, tv->tv+(alpha0dual*tv)*alpha0);
if TMTTmult(R0, GramS11)<>GramS11 then beep(99911); fi;
Add(sroots, alpha0);
Add(refls, R0);


for ii in [1..9] do 
  talpha:=List(zv);
  talpha[ii+1]:=1; talpha[ii+2]:=-1; 
  talphadual:=talpha*GramS11;
  if talphadual*talpha<>-2 then beep(88229191); fi;
  tR:=List(idmat11, tv->tv+(talphadual*tv)*talpha);
  if TMTTmult(tR, GramS11)<>GramS11 then beep(9119911); fi;
  Add(sroots, talpha);
  Add(refls, tR);
od;

talpha:=List(zv);
talpha[11]:=1;
talphadual:=talpha*GramS11;
if talphadual*talpha<>-1 then beep(55223391); fi;
tR:=List(idmat11, tv->tv+2*(talphadual*tv)*talpha);
if TMTTmult(tR, GramS11)<>GramS11 then beep(9119911); fi;
Add(sroots, talpha);
Add(refls, tR);


talpha:=List(zv);
talpha[1]:=3;
for kk in [2..11] do talpha[kk]:=-1;od;
talphadual:=talpha*GramS11;
if talphadual*talpha<>-1 then beep(22223391); fi;
tR:=List(idmat11, tv->tv+2*(talphadual*tv)*talpha);
if TMTTmult(tR, GramS11)<>GramS11 then beep(9119911); fi;
Add(sroots, talpha);
Add(refls, tR);


if Rank(sroots)<>11 then beep(56612); fi;
if CokerTorsion(sroots)<[] then beep(76612); fi;


GramL11:=2*GramS11;

discL11:=DiscriminantForm(GramL11);

vsdiscL11:=VsDiscg(discL11.discg);

qrefls:=List(refls, tg->OLtoOqL(tg, discL11));
qreflsperm:=List(qrefls, tqg->OqLtoOqLPerm(tqg, discL11.discg, vsdiscL11));

Printn(Size(Group(qreflsperm)));

gap> Read("Task20260905Ub.g");
46998591897600 

#####