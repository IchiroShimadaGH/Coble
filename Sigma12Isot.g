#Read("Sigma12Isot.g");

discrecS12:=DiscriminantForm(tGram);

nn := 13;

Df:=2*discrecS12.discf;


# 整数 m (0 <= m < 2^13) を F_2^9 のベクトルに変換

MaskToVector := function(m)
    local ii, v;
    if m>=2^nn then buzz(77162); fi;
    v:=List([0..nn-1],ii -> (QuoInt(m, 2^ii) mod 2));
    return (v);
end;

IsIsotropicMask := function(m)
  local vv;
  vv:=MaskToVector(m);
  return( (vv*Df*vv) mod 4=0);
end;

# 非零 isotropic vectors

IsoVectorsS12 :=Filtered([1..2^nn-1], IsIsotropicMask);
IsoVectorsS12:=Set(IsoVectorsS12);

#############################################################################
# 全 totally isotropic subspaces の列挙
# 部分空間は、その全要素を整数に encode した昇順リストで表す。
#############################################################################

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


EvenOverS12Codes := function()
  local levels, alllevels, newlevels, 
  S, v, T, d, Splusv, x;
  #
  levels := [[0]];
  alllevels := [levels];
  #
  for d in [1..nn] do
    newlevels := [];
    for S in levels do
      for v in IsoVectorsS12 do
        if not v in S then
          Splusv:=List(S, x -> MaskXor(x, v));
          T:=Union(S, Splusv);
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
    Printn(d, Length(newlevels));
    levels := newlevels;
  od;
  return alllevels;
end;

CodesByDimensionS12 := EvenOverS12Codes();

List(CodesByDimension, Length);

savedata(IsoVectorsS12);
savedata(CodesByDimensionS12);

########