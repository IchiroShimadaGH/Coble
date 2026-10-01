

#Read("Task20260921Ua.g");


# even_lattice_genus_README2.md により詳しい使い方が書いてある．
Read("even_lattice_genus.g");


#r := EvenLatticeGenus([2,17],[30],[[17/30]]);;
# EvenLatticeGenus(signature, discriminant group, discriminant form)
## 出力
# ```gap
# rec(
#     count := h,
#     grams := [G1,...,Gh],
#     method := "...",
#     spinorPrimes := [...]
# )
# ```

# `grams` は genus 内の各同型類からちょうど一つの Gram 行列を返します。
# 同型は通常の整格子の同型で、向きや判別群の marking は固定しません。
# Gram 行列は整数対称行列で対角が偶数です。二次多項式の係数行列ではなく、
# 双線形形式の Gram 行列です。

# 数学的に存在しない場合だけ `count=0; grams=[];` となります。
# ソフトウェアのエラー、未対応 API、資源不足を「0類」と扱うことはありません。
# 計算失敗時には入力一時ファイルを残して GAP のエラーにします。

# `method` は用いた方法、`spinorPrimes` は spinor 法で使った隣接素数です。
# `0` は最初の代表を意味し、素数0で隣接格子を作るわけではありません。
# 代表行列は正準形ではなく、実行や SageMath の版によって変わり得ます。


Read("SplitConsTools.g");

readdata("specialadjsf");

nopss:=List(specialadjsf, xx->xx[1]);;
Sort(nopss);;
nopss:=Reversed(nopss);;

Printn(Length(specialadjsf), List([1..10], ii->nopss[ii]));

cc:=0;
for stadj in specialadjsf do
  cc:=cc+1;
  tadj:=stadj[2];
  trec:=WGraphToGramh(tadj, Length(tadj));
  tGram:=trec.Gram;
  ttmdiscrec:=DiscriminantForm(-tGram);
  tsign:=[22, 3, 19]-SignatureQ(tGram);
  ttsign:=[tsign[2], tsign[3]];
  Trec:=EvenLatticeGenus(ttsign,ttmdiscrec.discg,ttmdiscrec.discf);;
  if Trec.count=0 then buzz(47652);
  else 
    Printn(cc, stadj[1], Length(tadj), Length(tGram),tsign[1],  ttmdiscrec.discg);
    if Trec.count>1 then Printn("h>1", Trec.count); fi;
  fi;
od;






############