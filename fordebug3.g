
readdata("ttbasisrec");
Read("t3.g");

t3GramL:=ttbasisrec.Gram;
t3basisrec:=BasisRec(t3GramL, 3);
NewOGLat(t3basisrec);