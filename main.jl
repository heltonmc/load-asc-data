include("load_convolve_data.jl")
include("DE-infinite-TR.jl")
using LsqFit
using Plots

IRF_filename = "IRF_10mm_SDS_750nm.asc"
DTOF_filename = "IL5_10mm_SDS_750nm_20ml.asc"
IRF = loadData(IRF_filename)
DTOF = loadData1(DTOF_filename)

function fitDTOF(DTOF,IRF)

    maxval,maxin = findmax(DTOF.counts)

    function conv_DT_IRF(t,β)
        #function we are fitting raw data to on line 30
        ρ=1
        Rt = (time = t, counts = refl_DT1(t,β,ρ))
        convDT = convolveTR_IRF(Rt,IRF)
        ~,maxinConv = findmax(convDT.counts)
        timefit = convDT.time .- convDT.time[maxinConv]
        log.(convDT.counts[maxinConv-maxin:maxinConv+(length(DTOF.counts)-maxin-1)])
    end

    lb = [0.001, 1]
    ub = [0.5, 80]
    w  = Vector(1:1:length(DTOF.counts)).^1
    β0 = [0.3, 40]
    fit = curve_fit(conv_DT_IRF, IRF.time[1:length(DTOF.time)], log.(DTOF.counts./maximum(DTOF.counts)),w,β0,lower = lb,upper = ub)
end

maxval,maxin = findmax(DTOF.counts)
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
