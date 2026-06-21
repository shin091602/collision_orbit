using DataFrames
using Arrow

# 1. 読み込みたいArrowファイルのパスを指定します
# （例：ヤコビ定数 2.94 の時系列データ）
file_path = joinpath(@__DIR__, "../../../data/CL_ver2/planar/timeseries/timeseries_C_2.941.arrow")
meta_file_path = joinpath(@__DIR__, "../../../data/CL_ver2/planar/metadata_all.arrow")
ana_file_path = joinpath(@__DIR__, "../../../data/CL_ver2/planar/metadata_analyzed.arrow")
neg_kep_file_path = joinpath(@__DIR__, "../../../data/CL_ver2/planar/negative_kepler_fractions.arrow")

# 2. データを読み込んでDataFrameに変換します
df_time = DataFrame(Arrow.Table(file_path))
df_meta = DataFrame(Arrow.Table(meta_file_path))
df_ana = DataFrame(Arrow.Table(ana_file_path))
df_neg_kep = DataFrame(Arrow.Table(neg_kep_file_path))
# 3. vscodedisplay でVS Codeのビューアを起動します
if @isdefined vscodedisplay
    # vscodedisplay(df_time)
    # vscodedisplay(df_meta)
    # vscodedisplay(df_ana)
    vscodedisplay(df_neg_kep)
else
    display(df_time)
    display(df_meta)
end