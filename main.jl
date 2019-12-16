include("load_asc_data.jl")
include("DE-infinite-TR.jl")

IRF = loadData("/Users/michaelhelton/Documents/load-asc-data/IRF_10mm_SDS_750nm.asc")
DTOF = loadData("/Users/michaelhelton/Documents/load-asc-data/IL5_10mm_SDS_750nm_c1.asc")
DTOF1 = loadData1("/Users/michaelhelton/Documents/load-asc-data/IL5_10mm_SDS_750nm_c1.asc")

maxval,maxin = findmax(DTOF1.counts)
time1 = []


function conv_DT_IRF(t,β)
    Rt = (time = t, counts = refl_DT(t,β))
    convDT = convolveTR_IRF(Rt,IRF)
    ~,maxindex = findmax(convDT.counts)
    timefit = convDT.time .- convDT.time[maxindex]
    convDT.counts[maxindex-maxin:maxindex+(length(DTOF1.counts)-maxin-1)]
end

using LsqFit
lb = [0.01, 1]
ub = [0.5, 40]
w  = Vector(1:1:length(DTOF1.counts)).^4

fit = curve_fit(conv_DT_IRF, IRF.time[1:length(DTOF1.time)], DTOF1.counts./maximum(DTOF1.counts),w,[0.1, 10],lower = lb,upper = ub)


beta_fit = fit.param

# preparing the fitting evaluation
xfit = IRF.time[1:length(DTOF1.time)]
yfit = conv_DT_IRF(xfit,fit.param)

# Creating a new figure object

# Plotting two datasets
plot(DTOF1.counts./maximum(DTOF1.counts), color="black", linewidth=2.0,yaxis =:log)
plot!(yfit, color="red", linewidth=2.0)

# Labeling axes
#xlabel("x", fontsize="xx-large")
#ylabel("y", fontsize="xx-large")
