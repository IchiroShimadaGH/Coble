#Read("Task20260910Uc.g");


readdata("Rrecs");
# ogsize の大きい順に並べる

DD:=46998591897600;

indexes:=List(Rrecs, xx->DD/xx.imsize);

SortParallel(indexes, Rrecs);

# TeX table の各行を出力
for numb in [1..Length(Rrecs)] do
    Printn(
        numb, " & ",
        ADElatex(Rrecs[numb].ADEtype), " & ",
        Rrecs[numb].ogsize, " & ",
        DD/Rrecs[numb].imsize, " & ",
        Rrecs[numb].ms[1], " & ",
        Rrecs[numb].ms[2],
        " \\\\"
    );
od;


####################