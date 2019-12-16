using DelimitedFiles

function loadData(filename)
    data_array = readdlm(filename,skipstart=10)
    t = data_array[1:end-1,1]
    Rt = data_array[1:end-1,2]
    (time = convert(Array{Float64},t), counts = convert(Array{Float64},Rt))
end

function process(data)
    time = data.time
    counts = data.counts

    maxvalue, maxindex = findmax(counts)
    ind = maxindex
    while counts[ind] > 0.01*maxvalue || counts[ind] > counts[ind-1]
        ind -= 1
    end

    noise = counts[ind]
    counts = counts .- noise
    inds = counts .> 10*noise
    (time = time[inds], counts = counts[inds])
end

function loadData1(filename)
    data_array = readdlm(filename,skipstart=10)
    time = data_array[1:end-1,1]
    counts = data_array[1:end-1,2]

    maxvalue, maxindex = findmax(counts)
    ind = maxindex
    while counts[ind] > 0.01*maxvalue || counts[ind] > counts[ind-1]
        ind -= 1
    end

    noise = counts[ind]
    counts = counts .- noise
    inds = counts .> 10*noise
    (time = (time[inds].-time[maxindex]), counts = counts[inds])
end

using DSP
using Statistics

function convolveTR_IRF(TR,IRF)

    convolved = conv(IRF.counts,TR.counts)
    convolved = convolved./maximum(convolved)
    #find new time vector

    n = collect(1:length(TR.time)-1)
    tbin = zeros(length(TR.time))
    for a in n
       tbin[a] = TR.time[a+1]-TR.time[a]
    end
    tavg = mean(tbin[1:end-1])
    (time = collect(0:1:length(convolved)-1).*tavg, counts = convolved)

end
