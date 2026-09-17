# 不定値偶格子の genus：存在・類数・全代表

## 必要なもの

- GAP 4
- SageMath（10.9 で検証済み）
- `even_lattice_genus.g` と `lattice_genus_backend.py` を同じディレクトリに置く

追加の GAP パッケージ、Magma、インターネット接続は不要です。
SageMath の一部の内部 API を使用しているため、他の版では付属テストを実行してください。

## 呼び出し

```gap
Read("even_lattice_genus.g");
r := EvenLatticeGenus([2,17],[30],[[17/30]]);;
r.count;     # 1
r.grams[1];  # 19 x 19 の偶整対称 Gram 行列
```

SageMath が PATH にない場合は、実行前に次を設定します。

```gap
EvenLatticeGenusSage := "/absolute/path/to/sage";
```

`Read` に GAP ファイルの絶対パスを指定しても動きます。
バックエンドの自動検出を上書きする場合は次を設定します。

```gap
EvenLatticeGenusBackend := "/absolute/path/to/lattice_genus_backend.py";
```

## 入力の規約

### SageMath から直接使う場合

GAP を使わず、同じディレクトリにある Python モジュールを読み込みます。
SageMath のコンソール／ノートブック／`.sage` ファイルで：

```python
from lattice_genus_backend import even_genus_representatives
grams = even_genus_representatives([2,17], [30], [[17/30]])
len(grams)  # 類数：1
grams[0]   # Sage の ZZ 上の行列
```

非存在の場合は空リスト `[]` が返り、類数は `len(grams) == 0` です。
判別形式にはネストしたリストでも `matrix(QQ, ...)` でも指定できます。
Sage の preparser を通さない通常の Python コードでは、`17/30` は浮動小数点数に
なるため、必ず `QQ(17)/30` のように正確な有理数を作ってください。
別ディレクトリのノートブックからは、モジュールのあるディレクトリを
`sys.path` に追加してから import してください。
`lattice_genus_examples.sage` にそのまま実行できる例があります。

### 引数

```gap
EvenLatticeGenus([splus,sminus],invariants,form)
```

- `splus,sminus` はともに正の整数。不定値・非退化の格子が対象です。
- `invariants = [d1,...,dm]` は、各生成元の位数を表すリストです。
  判別群は直和 `Z/d1 + ... + Z/dm`、各 `di >= 2` です。
  不変因子の整除順でなくても、直和分解として正しければ使用できます。
- `form` はこの生成元に関する対称有理数行列です。
  二次形式の値は **Q/2Z**、`q(x)=x*form*x` です。
  対角成分は modulo 2Z、非対角成分は modulo Z の情報を持ちます。
- 自明な判別群は `invariants=[]; form=[];` とします。
- 入力は近似小数でなく、`17/30` のような GAP の正確な有理数を使ってください。
- 判別形式の符号は入力どおり使います。自動的に反転しません。

生成元の位数と行列の整合性、二次形式の well-definedness、双線形形式の
非退化性も検証します。入力自体が不正ならエラーになります。

## 出力

```gap
rec(
    count := h,
    grams := [G1,...,Gh],
    method := "...",
    spinorPrimes := [...]
)
```

`grams` は genus 内の各同型類からちょうど一つの Gram 行列を返します。
同型は通常の整格子の同型で、向きや判別群の marking は固定しません。
Gram 行列は整数対称行列で対角が偶数です。二次多項式の係数行列ではなく、
双線形形式の Gram 行列です。

数学的に存在しない場合だけ `count=0; grams=[];` となります。
ソフトウェアのエラー、未対応 API、資源不足を「0類」と扱うことはありません。
計算失敗時には入力一時ファイルを残して GAP のエラーにします。

`method` は用いた方法、`spinorPrimes` は spinor 法で使った隣接素数です。
`0` は最初の代表を意味し、素数0で隣接格子を作るわけではありません。
代表行列は正準形ではなく、実行や SageMath の版によって変わり得ます。

## 例

以前の signature (2,15) の例：

```gap
r := EvenLatticeGenus(
    [2,15], [3,3,42],
    -[[4/3,0,2/3],[0,0,1/3],[2/3,1/3,25/42]]
);;
r.count;  # 1
```

非存在：

```gap
r := EvenLatticeGenus([1,2],[30],[[17/30]]);;
r.count;  # 0
r.grams;  # []
```

階数2で二つの同型類：

```gap
r := EvenLatticeGenus([1,1],[229],[[-2/229]]);;
r.count;  # 2
r.grams;
```

階数3で二つの同型類：

```gap
r := EvenLatticeGenus(
    [2,1], [2,34,578],
    [[1/2,0,0],[0,1/34,0],[0,0,-1/578]]
);;
r.count;  # 2
```

階数3で四つの同型類：

```gap
p := 1513;;
r := EvenLatticeGenus(
    [2,1], [2,2*p,2*p^2],
    [[1/2,0,0],[0,1/(2*p),0],[0,0,-1/(2*p^2)]]
);;
r.count;  # 4
```

## アルゴリズムと完全性

1. SageMath の有限二次形式から、指定 signature での存在を判定し、genus を構成。
2. `rank >= length(discriminant group)+2` なら Nikulin の一意性を使用。
3. 階数2では、判別式 `abs(det)` の二元二次形式の全被約候補を列挙し、
   genus で絞り込み、GL(2,Z) 同値で重複を除去。
   内容が1でない形式と平方判別式も含みます。列挙境界は整数平方根で計算します。
4. その他の階数3以上では、SageMath の局所 spinor kernel を用い、
   improper spinor class の有限商群の**全要素**を列挙。
   各非自明剰余類を表す良い奇素数 p を選び、基準格子の p-neighbor を構成します。
   不定値・階数3以上では improper spinor class と通常の整同型類が一致します。
5. 全出力の偶性、行列式、signature を含む genus の一致を検証。

SageMath 10.9 の `spinor_generators()` は部分群更新時に商群の全生成元を含めて
しまい、非自明な素数を一つ得たところで停止し得ます。
また、生成元だけの隣接格子では複数生成元の積に対応する剰余類を網羅しません。
このプログラムは `G.representatives()` に列挙を丸投げせず、有限商群の
全剰余類を明示的に処理し、商群の位数と代表数の一致を確認します。
この回避策は同梱ファイル内だけで行い、インストール済み SageMath は変更しません。

## 検証

同じディレクトリで：

```sh
sage -python test_lattice_genus_backend.py
```

以前の2例、非存在、U、U(2)、二元の2類、三元の2類・4類、
不正入力の拒否を確認しています。二元の場合は平方判別式を含む
判別式80以下の全 genus について SageMath の別の公開列挙関数と類数を照合しています。
数学的な完全性は上記の分類アルゴリズムによるもので、有限範囲の探索から推測していません。

## 参考

- [SageMath: 有限二次形式、is_genus、genus](https://doc.sagemath.org/html/en/reference/modules/sage/modules/torsion_quadratic_module.html)
- [SageMath: genus と spinor genus の実装](https://github.com/sagemath/sage/blob/develop/src/sage/quadratic_forms/genera/genus.py)
- [SageMath: 二元二次形式](https://doc.sagemath.org/html/en/reference/quadratic_forms/sage/quadratic_forms/binary_qf.html)

一般には行列式・階数・類数が大きいほど計算が重くなります。
探索回数の上限で打ち切って不完全なリストを返す処理は入れていません。
