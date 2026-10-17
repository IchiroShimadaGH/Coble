
 #Read("Task20261017Ub.g");

readdata("newsingK3configs");

#  OLrec:=rec(
#       Gram:=newGram, 
#       h:=newth, 
#       basisdual:=newbasisdual,
#       extgroupgens:=extgroupgens, 
#       extdeg:=nopsxtgroup
#     );


new2singK3configs:=[];

for trec in newsingK3configs do
  Printn(trec.dd, trec.GammaSize, trec.hperprec.AutShsize);
  trec.nonsingOLs:=NonSingOverLats(trec.GramS, trec.h, trec.hperprec.AutShgens);
  Printn(Length(trec.nonsingOLs));
  for olrec in trec.nonsingOLs do 
    olrec.sprats:=GetSpRats(olrec.Gram, olrec.h);
    Printn("_____", olrec.extdeg, olrec.sprats.nopss);
  od;
  #
  Add(new2singK3configs, trec);
  savedata(new2singK3configs);
  Printn("_____________");
od;


# gap> Read("Task20261017Ub.g");
# 360 4455 11520 
# 1 
# _____ 1 [ 0, 4455 ] 
# _____________ 
# 392 4263 5376 
# 1 
# _____ 1 [ 0, 4263 ] 
# _____________ 
# 432 4059 1728 
# 1 
# _____ 1 [ 0, 4059 ] 
# _____________ 
# 440 4020 480 
# 1 
# _____ 1 [ 0, 4020 ] 
# _____________ 
# 468 3888 576 
# 1 
# _____ 1 [ 0, 3888 ] 
# _____________ 
# 488 3813 96 
# 1 
# _____ 1 [ 0, 3813 ] 
# _____________ 
