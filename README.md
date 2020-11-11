# 5G_NFFFT
5G near field to far field transformation algorithm
11/11に別レポジトリに移動

## 仕様
グラフ表示にgnuplotが必要

変数名などのスペルミス多数

## 使い方
makeすれば自動でbinディレクトリにfield_transという実行ファイルがコンパイルされると思います。

field_trans \<near field data\> \<far field data\>

とすれば変換結果が出力されます。標準出力に多くのログが出力されるので適当なテキストファイルにリダイレクトしたほうがいいです。

## data
dataディレクトリに自分が用いた近傍界、遠方界のデータが入っています。
現時点では厳密にこの形のテキスト形式のファイルしか読むことはできません。
ファイル内のnear field とかfar fieldとかz = ... のような記述は間違っている（というより関係がない）ので注意してください。

ファイル名は以下のように書いてあります。

\<antenna name>_ <走査面>_ <測定距離[mm]>_ <測定点数> .txt

例えば

array_polar_15_32.txtはアレーアンテナの球面走査で測定距離は15mmで３２点サンプリングしたデータです。

