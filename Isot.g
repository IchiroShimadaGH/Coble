#Read("Tak20260918Ua.g");

#############################################################################
# L = <2> + <-2>^8
#############################################################################

nn := 9;
DD := DiagonalMat(Concatenation([2], List([1..8], i -> -2)));

# 整数 m (0 <= m < 2^9) を F_2^9 のベクトルに変換

MaskToVector := function(m)
    local ii, v;
    if m>=2^nn then buzz(77162); fi;
    v:=List([0..nn-1],ii -> (QuoInt(m, 2^ii) mod 2) * One(GF(2))
    return (v);
end;

# x_0 - x_1 - ... - x_8 mod 4 を計算
QNumeratorMod4 := function(m)
    local xx;
    xx := List([0..nn-1], i -> QuoInt(m, 2^i) mod 2);
    return (xx[1] - Sum(xx{[2..nn]})) mod 4;
end;

IsIsotropicMask := function(m)
    return(QNumeratorMod4(m) = 0);
end;

# 非零 isotropic vectors
IsoVectors :=Filtered([1..2^nn-1], IsIsotropicMask);

#############################################################################
# 全 totally isotropic subspaces の列挙
#
# 部分空間は、その全要素を整数に encode した昇順リストで表す。
#############################################################################

EvenOverlatticeCodes := function()
    local levels, alllevels, newlevels, S, v, T, d;

    # 零部分空間
    levels := [[0]];
    alllevels := [levels];

    for d in [1..nn] do
        newlevels := [];

        for S in levels do
            for v in IsoVectors do
                if not v in S then

                    # <S,v> = S union (v+S)
                    T := Set(Concatenation(
                        S,
                        List(S, x -> LogXor(x, v))
                    ));

                    if ForAll(T, IsIsotropicMask) then
                        AddSet(newlevels, T);
                    fi;
                fi;
            od;
        od;

        if Length(newlevels) = 0 then
            break;
        fi;

        Add(alllevels, newlevels);
        levels := newlevels;
    od;

    return alllevels;
end;

CodesByDimension := EvenOverlatticeCodes();

List(CodesByDimension, Length);