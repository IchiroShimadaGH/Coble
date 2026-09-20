# Read("WkTools.g");

# 2026/09/26 読み直し


MinRepWk := function(kk, w)
  local Ek, n, edgepos, r, i, j, tau, images, mask, s, cand, e, im1, im2, pos, a, best;
  #
  if not IsInt(kk) or kk < 2 then
    Error("kk must be an integer greater than or equal to 2");
  fi;
	#
  # List of edges in the prescribed order.
  Ek := [];
  for i in [1..kk] do
    for j in [i+1..kk] do
      Add(Ek, [i, j]);
    od;
  od;
	#
  n := Length(Ek);
	#
  if not IsList(w)
     or Length(w) <> n
     or not ForAll(w, a -> a in [0..4]) then
    Error("w must be a list of length Binomial(kk,2)",
        " with entries in [0..4]");
  fi;
	#
  # edgepos[i][j] is the position of {i,j} in Ek.
  edgepos := List([1..kk], i -> List([1..kk], j -> 0)); #initialize
	#
  for r in [1..n] do
    i := Ek[r][1];
    j := Ek[r][2];
    edgepos[i][j] := r;
    edgepos[j][i] := r;
  od;
	#
  best := ShallowCopy(w);
	#
  for tau in SymmetricGroup(kk) do
		#
    # images[i] = i^tau
    images := List([1..kk], i -> i^tau);
		#
    # Simultaneously changing all signs has no effect.
    # Hence we may impose s[1] = 1.
		#
    for mask in [0..2^(kk-1)-1] do
			#
      s := List([1..kk], i -> 1);
      for i in [2..kk] do
        if QuoInt(mask, 2^(i-2)) mod 2 = 1 then
          s[i] := -1;
        fi;
      od;
			#
      cand := [];
			#
      for e in Ek do
        i := e[1];
        j := e[2];
        im1 := images[i];
        im2 := images[j];
        pos := edgepos[im1][im2];
        a := w[pos];
        if s[i] * s[j] = -1 then
          a := 4-a;
        fi;
        Add(cand, a);
      od;
      # GAP compares lists lexicographically.
      if cand < best then
        best := cand;
      fi;
    od;
  od;
  return best;
end;

MinRepsListWk := function(kk)
  local Ek, n, edgepos, i, j, r,total, seen, reps, Sk,encode, 
	decode,code, w, tau, images, mask, s, cand, e,im1, im2, pos, a, candcode;
	#
  if not IsInt(kk) or kk < 2 then
    Error("kk must be an integer greater than or equal to 2");
  fi;
	#
  # Edges in the prescribed order.
  Ek := [];
  for i in [1..kk] do
    for j in [i+1..kk] do
      Add(Ek, [i, j]);
    od;
  od;
	#
  n := Length(Ek);
	#
  # edgepos[i][j] is the position of {i,j} in Ek.
  edgepos := List([1..kk], i -> List([1..kk], j -> 0));# #initial setting 
	#
  for r in [1..n] do
    i := Ek[r][1];
    j := Ek[r][2];
    edgepos[i][j] := r;
    edgepos[j][i] := r;
  od;
	#
  # Encode [a_1,...,a_n] as a base-5 integer.
  # This encoding preserves lexicographic order.
  encode := function(v)
    local c, x;

    c := 0;
    for x in v do
      c := 5*c+x;
    od;

    return c;
  end;
  # Inverse of encode.
  decode := function(c)
    local v, q, t;

    v := List([1..n], x -> 0);
    q := c;

    for t in [n, n-1 .. 1] do
      v[t] := q mod 5;
      q := QuoInt(q, 5);
    od;
    return v;
  end;

  total := 5^n;

  # Compact Boolean list indexed by code+1.
  seen := BlistList([1..total], []);
  reps := [];

  Sk := SymmetricGroup(kk);

  # The codes increase in the lexicographic order of W_k.
  for code in [0..total-1] do
    if not seen[code+1] then

      w := decode(code);

      # Since this is the first element of its orbit encountered,
      # it is the lexicographically minimal representative.
      Add(reps, w);

      # Mark the entire G-orbit of w.
      for tau in Sk do
        images := List([1..kk], i -> i^tau);

        # We may fix s_1=1, since simultaneous reversal of
        # all signs acts trivially.
        for mask in [0..2^(kk-1)-1] do

          s := List([1..kk], i -> 1);

          for i in [2..kk] do
            if QuoInt(mask, 2^(i-2)) mod 2 = 1 then
              s[i] := -1;
            fi;
          od;

          cand := [];

          for e in Ek do
            i := e[1];
            j := e[2];

            im1 := images[i];
            im2 := images[j];

            pos := edgepos[im1][im2];
            a := w[pos];

            if s[i]*s[j] = -1 then
              a := 4-a;
            fi;

            Add(cand, a);
          od;

          candcode := encode(cand);
          seen[candcode+1] := true;
        od;
      od;
    fi;
  od;

  return reps;
end;

OrbitWk := function(kk, w)
  local Ek, n, edgepos, i, j, r,
      orbit, tau, images, mask, s,
      cand, e, im1, im2, pos, a;

  if not IsInt(kk) or kk < 2 then
    Error("kk must be an integer greater than or equal to 2");
  fi;

  # Edges in the prescribed order.
  Ek := [];
  for i in [1..kk] do
    for j in [i+1..kk] do
      Add(Ek, [i, j]);
    od;
  od;

  n := Length(Ek);

  if not IsList(w)
     or Length(w) <> n
     or not ForAll(w, a -> a in [0..4]) then
    Error("w must be a list of length Binomial(kk,2)",
        " with entries in [0..4]");
  fi;

  # edgepos[i][j] is the position of {i,j} in Ek.
  edgepos := List([1..kk],
          i -> List([1..kk], j -> 0));

  for r in [1..n] do
    i := Ek[r][1];
    j := Ek[r][2];
    edgepos[i][j] := r;
    edgepos[j][i] := r;
  od;

  orbit := [];

  for tau in SymmetricGroup(kk) do
    images := List([1..kk], i -> i^tau);

    # Fix s[1]=1, since simultaneous reversal of all signs
    # acts trivially.
    for mask in [0..2^(kk-1)-1] do

      s := List([1..kk], i -> 1);

      for i in [2..kk] do
        if QuoInt(mask, 2^(i-2)) mod 2 = 1 then
          s[i] := -1;
        fi;
      od;

      cand := [];

      for e in Ek do
        i := e[1];
        j := e[2];

        im1 := images[i];
        im2 := images[j];

        pos := edgepos[im1][im2];
        a := w[pos];

        if s[i]*s[j] = -1 then
          a := 4-a;
        fi;

        Add(cand, a);
      od;

      Add(orbit, cand);
    od;
  od;

  # Remove duplicates and sort lexicographically.
  return Set(orbit);
end;

#####