#Read("Task20260918Ua.g");

# Made by JeePeeChee

#############################################################################
# L = <2> + <-2>^8
#############################################################################

nn := 9;
DD := DiagonalMat(Concatenation([2], List([1..8], i -> -2)));
DDdual:=InverseMat(DD);

# 整数 m (0 <= m < 2^9) を F_2^9 のベクトルに変換

MaskToVector := function(m)
    local ii, v;
    #if m>=2^nn then buzz(77162); fi;
    v:=List([0..nn-1],ii -> (QuoInt(m, 2^ii) mod 2) * One(GF(2)));
    return (v);
end;

MaskXor := function(a, b)
    local ans, power, i;

    ans := 0;
    power := 1;

    for i in [1..nn] do
        if (a mod 2) <> (b mod 2) then
            ans := ans + power;
        fi;

        a := QuoInt(a, 2);
        b := QuoInt(b, 2);
        power := 2 * power;
    od;

    return ans;
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
                        List(S, x -> MaskXor(x, v))
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

Printn(List(CodesByDimension, Length));

###################

CodeBasisVectors := function(S)
    local VV;
    VV := VectorSpace(GF(2), List(S, MaskToVector));
    return BasisVectors(Basis(VV));
end;

OverlatticeGramMatrix := function(S)
    local CB, A, H, latrec;

    CB := CodeBasisVectors(S);

    # M は Z^9 と、codeword/2 で生成される．
    # 2M の生成行列を作る．
    A := Concatenation(
        2 * IdentityMat(nn),
        List(CB, v -> List(v, IntFFE))
    );

    H := HermiteNormalFormIntegerMat(A);
    H := Filtered(H, r -> not ForAll(r, IsZero));
    tGram:=(1/4)*TMTTmult(H, DD);
    if not IsIntMat(tGram) then buzz(51112); fi;
    latrec:=rec(
        basis:=H/2,
        basisdual:=H/2*DD, 
        Gram:=tGram:
    );

    return latrec;
end;

