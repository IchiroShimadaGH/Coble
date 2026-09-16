# Keep this file and lattice_genus_backend.py in the same directory.
# SageMath 10.9 tested; no optional GAP package required.
# q(x) = x*form*x mod 2Z. Ordinary isometry, not oriented/proper isometry.

EvenLatticeGenusBackend := (function()
    local name, positions;
    name := INPUT_FILENAME();
    positions := Positions(name,'/');
    if Length(positions)=0 then
        return Filename(DirectoryCurrent(),"lattice_genus_backend.py");
    fi;
    return Concatenation(name{[1..Last(positions)]},"lattice_genus_backend.py");
end)();
EvenLatticeGenusSage := Filename(DirectoriesSystemPrograms(),"sage");

EvenLatticeGenus := function(signature,invariants,form)
    local m, directory, request, output, stream, input, log, status, result;
    if not IsList(signature) or Length(signature)<>2
       or not ForAll(signature,x -> IsInt(x) and x>0) then
        Error("signature must be [splus,sminus], both positive");
    fi;
    if not IsList(invariants) or not ForAll(invariants,x -> IsInt(x) and x>=2) then
        Error("invariants must list cyclic orders >=2 (or [] for the trivial group)");
    fi;
    m := Length(invariants);
    if not IsList(form) or Length(form)<>m
       or not ForAll(form,row -> IsList(row) and Length(row)=m and ForAll(row,IsRat)) then
        Error("form must be a rational square matrix matching invariants");
    fi;
    if EvenLatticeGenusSage=fail or not IsExecutableFile(EvenLatticeGenusSage) then
        Error("Set EvenLatticeGenusSage to the absolute path of your sage executable");
    fi;
    if not IsReadableFile(EvenLatticeGenusBackend) then
        Error("Cannot read lattice_genus_backend.py; set EvenLatticeGenusBackend to its path");
    fi;
    directory := DirectoryTemporary();
    request := Filename(directory,"input.json");
    output := Filename(directory,"output.g");
    stream := OutputTextFile(request,false);
    SetPrintFormattingStatus(stream,false);
    PrintTo(stream,"{\"signature\":",String(signature),
        ",\"invariants\":",String(invariants),
        ",\"form\":",String(List(form,row -> List(row,String))),"}");
    CloseStream(stream);
    log := "";
    stream := OutputTextString(log,true);
    input := InputTextNone();
    status := Process(DirectoryCurrent(),EvenLatticeGenusSage,input,stream,
        ["-python",EvenLatticeGenusBackend,request,output]);
    CloseStream(input);
    CloseStream(stream);
    if status<>0 or not IsReadableFile(output) then
        Error("Sage computation failed (this is NOT a zero class number). Output: ",log,
              " Input retained at ",request);
    fi;
    result := ReadAsFunction(output)();
    if not IsRecord(result) or not IsBound(result.count) or not IsBound(result.grams)
       or result.count<>Length(result.grams) then
        Error("Invalid backend result");
    fi;
    RemoveFile(request);
    RemoveFile(output);
    return result;
end;
