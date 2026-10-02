# rank 19 のランダム偶格子探索（GAP 4）

追加パッケージ不要。ノルムは (v,v)=v G v^t（その半分ではない）。

## 実行

このディレクトリをカレントディレクトリにして GAP を起動し、次を実行する。

```gap
Read("random_rank19.g");;
R19Search(rec(target := 10000, file := "my_lattices.g"));
```

`target` は今回追加する格子数。`maxTrials` の既定値は 100000。
保存は発見ごとに追記する。既存ファイルを上書きしない。

```gap
# 別の乱数種でさらに 10000 個追加
R19Search(rec(target := 10000, seed := 12345,
              file := "my_lattices.g", resume := true));

# より広い形の Gram 行列を探索（Smith 標準形で巡回性を検査）
R19Search(rec(mode := "dense", target := 10000,
              maxTrials := 100000, seed := 98765,
              offBound := 2, diagonalMax := 20,
              file := "dense_lattices.g"));

# 結果を読み込む
Read("my_lattices.g");;
Length(R19Saved);
R19Saved[1].gram;
R19Saved[1].determinant;
R19Saved[1].discriminantInvariantFactors;
```

`resume` に同じ seed を使うと同じ候補を再生成する。違う候補が欲しい場合は seed を変える。
保存ファイルは GAP コードなので、自分で生成したファイルを読み込むこと。
同じ保存ファイルを複数のプロセスから同時に書かない。

## 三条件を満たす理由

### 共通

G は 19×19 の対称整数行列で、対角成分は偶数。
各行で `G[i][i] - Sum(j != i, AbsInt(G[i][j])) >= 4` を満たす。
不等式 `2|x_i x_j| <= x_i^2+x_j^2` より

    x G x^t >= Sum_i (G_ii - Sum_{j != i}|G_ij|) x_i^2
             >= 4 Sum_i x_i^2.

したがって正定値、rank 19、偶格子、非零ベクトルのノルムは 4 以上。
ノルム 2 がないことは短いベクトルの近似探索ではなく、この不等式で保証する。
`normLowerBound := 4` は下界であり、最小ノルムの実測値ではない。

### chain（標準、高速）

対角成分を {6,8,...,40} から独立に選び、隣接非対角成分を 1、他を 0 とする。
`diagonalMax` で対角成分の上限を変更できる（奇数なら直下の偶数）。
行 1,...,18、列 2,...,19 の 18 次小行列は、対角成分が 1 の下三角行列なので行列式 1。
従って 18 次 determinantal divisor は 1、Smith 不変因子は

    1,...,1, det(G).

よって L^vee/L は Z/det(G)Z。chain では毎回の Smith 計算を省略できる。
`R19Cyclic(G)` で独立に検査することもできる。

### dense

非対角成分を [-offBound..offBound] から独立に選ぶ。
対角成分は行の絶対値和に 4 を足した値以上の最小偶数に、ランダムな非負偶数を追加する。
この追加分は 0,...,2*floor((diagonalMax-6)/2)。
**dense では diagonalMax は対角成分そのものの上限ではない。**
SmithNormalFormIntegerMat により最初の 18 不変因子がすべて 1 の候補だけ採用する。

## 重複と探索範囲

標準では重複除去をせず、同型な格子や同一の Gram 行列も保存する。
`distinctDet := true` を指定すると、一つの保存ファイルで行列式が異なる候補だけを追加する。
この設定では格子は互いに非同型になるが、同じ行列式を持つ別の格子も捨てる。
異なる保存ファイルの間では重複排除しない。
進捗表示の `total` は保存件数。

この探索分布は全格子上の一様分布ではなく、小さい判別式や最小ノルムがちょうど 4 の格子を優先しない。
条件を満たす例を大量に集める用途に向く。
固定したパラメータでは候補集合は有限なので、頭打ちになったら diagonalMax を増やす。
そもそも条件を満たす非同型格子は無限に存在する。例えば chain の対角をすべて 6 とし、
最初の対角だけ 6+2t にすると、行列式は t に対して正の傾きで増加する。

## 動作確認

GAP 上で chain 100 個と dense 30 個を生成し、全例の Smith 標準形を検査した。
行列式と dense の対角条件、および既存ファイルへの追加も確認した。
`sample_rank19.g` は chain 方式で生成した 1000 個の例。
