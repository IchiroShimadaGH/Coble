Read("IsNewSh.g");


GetGramP:=function(GramS, h)
  local hdual, GramP;
  hdual:=h*GramS;
  GramP:=TransposedMat([hdual])*[hdual]-GramS;
  if h*GramP*h<>2 then beep(998912); fi;
  return();
end;

GetIniData:=function(GramS, h)
  GramP:=GetGramP
  inirec:=rec(
    rank:=
    GramS:=GramS,
    GramP:=GramP,
    h:=h,

  )
  return(inirec);
end;

######