for grec in gammarecs do 
  n:=grec.n;
  if n=24 then continue; fi;
  #
  if n=8 then erec:=erecs[1];minellrec:=minellrecs[1];
  elif n=12 then erec:=erecs[2];minellrec:=minellrecs[2];
  elif n=16 then erec:=erecs[3];minellrec:=minellrecs[3];
  else beep(771626);
  fi;
  #
  vep:=grec.epsilon;
  tau:=grec.tau;
  Xitau:=grec.plusrec.Xitau;
  discmrec:=minellrec.discmrec;
  #
  nn:=function(xi)
    local tilxi;
    tilxi:=xi-tau/2;
    return(tilxi*GramLeech*tilxi);
  end;
  #
  #
  for xi in Xitau do 
    nnxi:=nn(xi);
    for xxi  in Xitau do 
      nnxxi:=nn(xxi); 
      delta:=SolutionMat(erec.minusrec.basis, xxi-xi)*discmrec.proj;
      tv:=ModDiscg(delta, discmrec.discg);
      tminell:=minellrec.minells[SinglePosition(minellrec.vs, tv)];
      flagA:=(nnxi<>nnxxi);
      flagB:=(tminell<>0);
      flagC:=(tminell+nnxi+nnxxi >=4);
      if (flagA or flagB) and (not flagC) then beep(999911); fi;
    od;
  od;
  Printn("done", n);
od;

