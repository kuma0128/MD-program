# MD-program

分子動力学・Monte Carlo計算で使う基本アルゴリズムを、C++17で実装した小さなライブラリです。
旧Fortranコードが参照していた未収録のグローバルモジュールを廃止し、入力値と状態を明示的に
渡せるAPIへ置き換えています。

## 実装内容

- `EWALD`: 周期境界条件下の実空間・逆空間Ewald和、自己エネルギー、Lennard-Jones項と力
- `MC`: Lennard-Jones粒子のNVT Metropolis Monte CarloとRDF集計
- `SHAKE`: velocity Verletに対する距離拘束
- `RATTLE`: 位置・速度の距離拘束を含む2段階積分
- `ROTATE`: Z-Y-Z Euler回転、quaternion回転、一様ランダム回転
- `GAUSS`: 極形式Box-Muller法による標準正規乱数とヒストグラム生成CLI

ベクトル、周期境界条件、拘束条件などの共通型は `MD/common.hpp` にあります。
すべての実数計算には `double` を使用します。

## ビルドとテスト

Makeを使う場合:

```sh
make
make test
```

CMakeを使う場合:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

テストはGaussian乱数の統計量、回転の距離保存、Ewaldの力とエネルギー勾配、Monte Carloの
エネルギー整合性、RDF、SHAKE/RATTLEの拘束条件を検証します。

## Gaussian乱数の実行例

```sh
./build/gaussian --samples 1000000 --bins 200 --width 0.02 --seed 80 \
  --output output.d
```

出力は各行が「bin中央値、サンプル数」で、最終行に平均と分散を記録します。オプションを
省略した場合は上記と同じ既定値を使います。

## API使用例

```cpp
#include "MC/monte_carlo.hpp"

std::vector<md::Vec3> positions{{1, 1, 1}, {3, 1, 1}};
md::MonteCarloParameters parameters;
parameters.box_length = 10.0;
parameters.cutoff = 4.0;
parameters.beta = 1.0;
parameters.max_displacement = 0.1;

md::MonteCarloNVT simulation(positions, parameters, 80);
simulation.sweep();
```

各APIは不正な寸法・質量・box長などを例外で報告します。EwaldとMonte Carloのcutoffは
minimum image conventionが成立するよう、box長の半分以下にしてください。
