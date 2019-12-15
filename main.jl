include("load_asc_data.jl")
include("DE-infinite-TR.jl")

IRF = loadData("/Users/michaelhelton/Documents/load-asc-data/IRF_10mm_SDS_750nm.asc")
DTOF = loadData("/Users/michaelhelton/Documents/load-asc-data/SolidBlock_10mm_SDS_750nm.asc")




function conv_DT_IRF(t,β)
    Rt = (time = t, counts = refl_DT(t,β))
    convDT = convolveTR_IRF(Rt,IRF)
    convDT.counts[1:length(DTOF.counts)]
    convDT.counts[1400:2000]
end

using LsqFit
lb = [0.01, 1]
ub = [0.5, 15]
fit = curve_fit(conv_DT_IRF, IRF.time[1400:2000], DTOF.counts[1400:2000], [0.1, 10],lower = lb,upper = ub)
