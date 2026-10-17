 #Read("Task20261017Ua.g");


#Set(singK3configs, RecNames)
#Read("NonSingOverLats.g");
#Read("MinimalWgs.g");

#readdata("singK3configs");

#the items of singK3configs[1]
#[ "h", "GramS", "rr", "sprats", "tt", "dd", "GammaSize", "LambdaSize", "gtypes" ]

# the output of HperpRec
#  trec:=rec(
#     rank:=Length(h),
#     GramS:=GramS,
#     h:=h,
#     GramN:=GramN,
#     T:=T,
#     Tinv:=Tinv,
#     Nbasisrec:=Nbasisrec,
#     AutShsize:=Nogrec.size, 
#     AutShgens:=OGShgens
#   );


newsingK3configs:=[];


for trec in singK3configs do
  Printn(trec.dd, trec.GammaSize);
  trec.hperprec:=HperpRec(trec.GramS, trec.h);
  Add(newsingK3configs, trec);
  Printn(trec.hperprec.AutShsize);
  savedata(newsingK3configs);
  Printn("______________");
od;





# gap> Read("Task20261017Ua.g");
# 360 4455 
# 11520 
# ______________ 
# 392 4263 
# 5376 
# ______________ 
# 432 4059 
# 1728 
# ______________ 
# 440 4020 
# 480 
# ______________ 
# 468 3888 
# 576 
# ______________ 
# 488 3813 
# 96 
# ______________ 




 ##############