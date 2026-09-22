

#Read("Task20260922Ua.g");

readdata("Sigma12rec");

GramS12:=Sigma12rec.Gram;


if DeterminantIntMat(GramS12)<>1 then beep(615261); fi;
svs:=ShortestVectors(GramS12, 2);;

if Collected(svs.norms)<>[ [ 2, 132 ] ] then buzz(51612); fi;


Read("SplitConsTools.g");

t0h:=MakeVectei(13, 1);
t0Gram:=DiagonalMats([ [[2]], (-2)*GramS12]);

Printn(IsGeom(t0Gram, t0h));






Read("even_lattice_genus.g");


beep(888811111); 

uptokk:=10;
savename:="result20260921V";

theresult:=[];

thetask:=function(kk, wg)
  local addvs, addv, newwg, trec, tGram, th, isgeom, minflag, scons,
  nooverlatflag;
  addvs:=Getaddvs(kk);
  for addv in addvs do 
    newwg:=CopyAppend(wg, addv);
    trec:=NewWGraphToGramh(newwg, kk+1);
    tGram:=trec.Gram;
    th:=trec.h;
    isgeom:=IsGeom(tGram, th);
    if isgeom=true then 
      minflag:=IsMinimal(kk+1, newwg);
      if minflag<>false then 
        scons:=GetTotalSplcons(tGram, th);
        nooverlatflag:=NoOverlattice(tGram, th);
        trec.wg:=newwg;
        trec.scons:=scons;
        trec.stabssize:=Length(minflag[2]);
        trec.nooverlatflag:=nooverlatflag;
        Add(theresult, trec);
        Printn("___no", Length(theresult));
        Printn("kk", kk+1, "rho", Length(tGram),
               "wgraph", Collected(newwg), "scons", Length(scons));
        savedataas(theresult, savename);
        if kk+1<=uptokk then
          thetask(kk+1, newwg);
        fi;
      fi;
    fi;
  od;
end;

uptokk:=5;
savename:="resultrho7";
thetask(1, []);


############