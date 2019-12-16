include("load_convolve_data.jl")
include("DE-infinite-TR.jl")

IRF = loadData("/Users/michaelhelton/Documents/load-asc-data/IRF_10mm_SDS_750nm.asc")
DTOF = loadData("/Users/michaelhelton/Documents/load-asc-data/IL5_10mm_SDS_750nm_c1.asc")
DTOF1 = loadData1("/Users/michaelhelton/Documents/load-asc-data/IL5_10mm_SDS_750nm_20ml.asc")

maxval,maxin = findmax(DTOF1.counts)

using LsqFit
lb = [0.001, 3]
ub = [0.5, 13]
w  = Vector(1:1:length(DTOF1.counts)).^4
β0 = [0.2, 10]
fit = curve_fit(conv_DT_IRF, IRF.time[1:length(DTOF1.time)], DTOF1.counts./maximum(DTOF1.counts),w,β0,lower = lb,upper = ub)


xfit = IRF.time[1:length(DTOF1.time)]
yfit = conv_DT_IRF(xfit,fit.param)

scatter(DTOF1.counts./maximum(DTOF1.counts), color="black",yaxis =:log,label = "Expt",markersize=3,alpha=0.8)
plot!(yfit, color="red", linewidth=3, label = "DT-fit",alpha = 0.8)
xlabel!("time (ns)")
ylabel!("Reflectance [counts]")
title!("\\mu_{a} = $(round(fit.param[1],digits=4)) cm^{-1}, \\mu_{s}' = $(round(fit.param[2],digits=3)) cm^{-1}")
#savefig("fit_plot.png")
#savefig("fit_plot.pdf")
