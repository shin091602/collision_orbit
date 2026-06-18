using DataFrames
using Arrow

# 1. 読み込みたいArrowファイルのパスを指定します
# （例：ヤコビ定数 2.94 の時系列データ）
file_path = joinpath(@__DIR__, "../../../data/CL_ver2/planar/timeseries/timeseries_C_2.94.arrow")

# 2. データを読み込んでDataFrameに変換します
df_time = DataFrame(Arrow.Table(file_path))

# 3. vscodedisplay でVS Codeのビューアを起動します
if @isdefined vscodedisplay
    vscodedisplay(df_time)
else
    display(df_time)
end