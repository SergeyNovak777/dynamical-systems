function get_local_max(timeseries)
    values, times = timeseries;
    values_len = length(values);
    local_maxs = Float64[];
    times_local_maxs = Float64[];
    for index in range(2, values_len-1)
        if values[index-1] <= values[index] >= values[index+1]
            push!(local_maxs, values[index]);
            push!(times_local_maxs, times[index]);
        end
    end
    return local_maxs, times_local_maxs
end

function get_local_min(timeseries)
    values, times = timeseries;
    values_len = length(values);
    local_mins = Float64[];
    times_local_mins = Float64[];
    for index in range(2, values_len-1)
        if values[index-1] >= values[index] <= values[index+1]
            push!(local_mins, values[index]);
            push!(times_local_mins, times[index]);
        end
    end
    return local_mins, times_local_mins
end

function drop_artifacts(array_local_maxs, array_local_mins)
    local_maxs, times_local_maxs = array_local_maxs;
    local_mins, times_local_mins = array_local_mins;

    if t_local_maxs[1] < t_local_mins[1]
        popfirst!(local_maxs)
        popfirst!(times_local_maxs)
    end
    if t_local_maxs[end] > t_local_mins[end]
        pop!(local_maxs)
        pop!(times_local_maxs)
    end
end

function check_spikes(array_local_maxs, array_local_mins)
    local_maxs, times_local_maxs = array_local_maxs;
    local_mins, times_local_mins = array_local_mins;
    
    len_local_mins = length(times_local_mins);
    index_local_mins = 1;
    index_local_maxs = 1;

    while index_local_mins != len_local_mins
        if times_local_maxs
    end
end