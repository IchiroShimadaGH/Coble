#Read("ConcatrandomPro.g");

IsGkMinimal:=function(kk, wg, absstab)
  #
  local tSwgens, minflag, Swtask, ttau, newwg;
  #
  tSwgens:=GetEkkRec(kk).Swgens;
  #
  minflag:=true;
  Swtask:=function(twg, jj)
    local newtwg, tpos;
    if jj>kk then return(); fi;
    #
    Swtask(twg, jj+1);
    if minflag=false then return(); fi;
    #
    newtwg:=List(twg);
    for tpos in tSwgens[jj] do 
      newtwg[tpos]:=-newtwg[tpos];
    od;
    if newtwg<wg then minflag:=false; return(); fi;
    if newtwg<>twg then 
      Swtask(newtwg, jj+1);
      if minflag=false then return(); fi;
    fi;
    #
  end;
  #
  for ttau in absstab do 
    newwg:=ActionTauWg(ttau, wg);
    if newwg<wg then minflag:=false; break; fi;
    Swtask(List(newwg), 2);
    if not minflag then break; fi;
  od;
  return(minflag);
end;