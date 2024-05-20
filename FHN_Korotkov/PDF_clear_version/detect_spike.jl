function get_local_max(data)
    x_range, t_range = data
    array_local_max = Float64[]
    array_t_local_max = Float64[]
    for index in range(2, length(x_range)-1 )
        if x_range[index-1] < x_range[index] > x_range[index+1]
            push!(array_local_max, x_range[index])
            push!(array_t_local_max, t_range[index])
        end
    end
    return array_local_max, array_t_local_max
end

function get_local_min(data)
    x_range, t_range = data
    array_local_min = Float64[]
    array_t_local_min = Float64[]
    for index in range(2, length(x_range)-1 )
        if x_range[index-1] > x_range[index] < x_range[index+1]
            push!(array_local_min, x_range[index])
            push!(array_t_local_min, t_range[index])
        end
    end
    return array_local_min, array_t_local_min
end

function drop_artifacts(data_local_max, data_local_min)
    local_maxs, t_local_maxs = data_local_max
    local_mins, t_local_mins = data_local_min

    if t_local_maxs[1] < t_local_mins[1]
        popfirst!(local_maxs)
        popfirst!(t_local_maxs)
    elseif t_local_maxs[end] > t_local_mins[end]
        pop!(local_maxs)
        pop!(t_local_maxs)
    end
end

function get_amplitude_spike(local_min_left, local_min_right, local_max)
    maxmin = maximum([local_min_left, local_min_right])
    amplitude = abs( local_max - maxmin )
    return amplitude
end

function get_amplitudes_all_events(array_local_max, array_local_min)
    amplitudes = Float64[]
    for index in range(1, length(array_local_max)-1)
        amplitude = get_amplitude_spike(array_local_min[index], array_local_min[index+1], array_local_max[index])
        push!(amplitudes, amplitude)
    end
    return amplitudes
end

function select_spikes(array_local_min, data_local_max, threshold)
    array_local_max, array_t_local_max = data_local_max
    array_local_max_above_thr = Float64[]
    array_t_local_max_above_thr = Float64[]
    amplitudes_above_thr = Float64[]
    for index in range(1, length(array_local_max)-1)
        amplitude = get_amplitude_spike(array_local_min[index], array_local_min[index+1], array_local_max[index])
        if amplitude >= threshold
            push!(array_local_max_above_thr, array_local_max[index])
            push!(array_t_local_max_above_thr, array_t_local_max[index])
            push!(amplitudes_above_thr, amplitude)
        end
    end
    return array_local_max_above_thr, array_t_local_max_above_thr, amplitudes_above_thr
end