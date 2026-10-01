#Read("Task20260917Ua.g");

Read("Wk.g");


GorbsW2:=MinRepsListWk(2);
Printn("k=2", Length(GorbsW2));
GorbsW3:=MinRepsListWk(3);
Printn("k=3", Length(GorbsW3));
GorbsW4:=MinRepsListWk(4);
Printn("k=4", Length(GorbsW4));
GorbsW5:=MinRepsListWk(5);
Printn("k=5", Length(GorbsW5));

savedata(GorbsW2);
savedata(GorbsW3);
savedata(GorbsW4);
savedata(GorbsW5);

# gap> Read("Task20260917Ua.g");
# k=2 3 
# k=3 14 
# k=4 161 
# k=5 6595


#####