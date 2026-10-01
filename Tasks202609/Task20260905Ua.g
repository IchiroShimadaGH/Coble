#Read("Task20260905Ua.g");


GramRdual:=InverseMat(GramR);
IntGramRdual:=2*GramRdual;
svsrec:=ShortestVectors(-IntGramRdual, 4);
Printn(Collected(svsrec.norms));
svsdual:=svsrec.vectors;
svs:=svsdual*GramRdual;


VsDiscg:=function(discg)
  local ii;
  return(Cartesian(List(discg, ii->[0..ii-1])));
end;

vsdiscR:=VsDiscg(discR.discg);




#####