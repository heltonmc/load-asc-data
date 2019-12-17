include("load_convolve_data.jl")
include("DE-infinite-TR.jl")
using LsqFit
using Plots

IRF_filename = "/Users/michaelhelton/Documents/load-asc-data/IRF_10mm_SDS_750nm.asc"
DTOF_filename = "/Users/michaelhelton/Documents/load-asc-data/IL5_10mm_SDS_750nm_20ml.asc"
IRF = loadData(IRF_filename)
DTOF = loadData1(DTOF_filename)
maxval,maxin = findmax(DTOF.counts)


function fitDTOF(DTOF,IRF)

    function conv_DT_IRF(t,β)
        ρ=1
        Rt = (time = t, counts = refl_DT1(t,β,ρ))
        convDT = convolveTR_IRF(Rt,IRF)
        ~,maxindex = findmax(convDT.counts)
        timefit = convDT.time .- convDT.time[maxindex]
        convDT.counts[maxindex-maxin:maxindex+(length(DTOF.counts)-maxin-1)]
    end

    lb = [0.001, 3]
    ub = [0.5, 50]
    w  = Vector(1:1:length(DTOF.counts)).^2
    β0 = [0.2, 10]
    fit = curve_fit(conv_DT_IRF, IRF.time[1:length(DTOF.time)], DTOF.counts./maximum(DTOF.counts),w,β0,lower = lb,upper = ub)
end

fit = fitDTOF(DTOF,IRF)
xfit = IRF.time[1:length(DTOF.time)]
yfit = conv_DT_IRF(xfit,fit.param)

scatter(xfit,DTOF.counts./maximum(DTOF.counts), color="black",yaxis =:log,label = "Expt",markersize=3,alpha=0.8)
plot!(xfit,yfit, color="red", linewidth=3, label = "DT-fit",alpha = 0.8)
xlabel!("time (ns)")
ylabel!("Reflectance [counts]")
title!("\\mu_{a} = $(round(fit.param[1],digits=4)) cm^{-1}, \\mu_{s}' = $(round(fit.param[2],digits=3)) cm^{-1}")
#savefig("fit_plot.png")
#savefig("fit_plot.pdf")
