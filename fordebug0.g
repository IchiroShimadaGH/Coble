#
PGUSTOS3:=[];
for tgS in PGUSOS3 do 
  Append(PGUSTOS3, List(PGUTOS3, tgT->tgS*tgT));
od;

if Length(PGUSTOS3)<>81*112 then beep(585871); fi;

PGUSTOS3g10:=List(PGUSTOS3, tg->tg*gdpp10);
PGUSTOS3g31:=List(PGUSTOS3, tg->tg*gdpp31);

####