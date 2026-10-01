#Read("ConcatrandomPro.g");


GetGramP:=function(GramS, h)
  local hdual, GramP;
  hdual:=h*GramS;
  GramP:=TransposedMat([hdual])*[hdual]-GramS;
  if h*GramP*h<>2 then beep(998912); fi;
  return(GramP);
end;



#############