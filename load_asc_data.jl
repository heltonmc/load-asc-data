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
    (time = time[inds], counts = counts[inds])
end
