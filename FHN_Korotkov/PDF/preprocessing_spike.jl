function check_timeseries(array_spikes_max, array_spikes_thresholds)
    if length(array_spikes_max)+1 == length(array_spikes_thresholds)
        flag = true
    else
        flag = false
    end
    return flag
end

function drop_false_start_end(array_t_spikes_max, array_spikes_max, array_t_spikes_thresholds)
    if array_t_spikes_max[1] < array_t_spikes_thresholds[1]
        popfirst!(array_t_spikes_max)
        popfirst!(array_spikes_max)
    elseif array_t_spikes_max[end] > array_t_spikes_thresholds[end]
        pop!(array_t_spikes_max)
        pop!(array_spikes_max)
    end
end

function get_amplitude_spike(local_min_left, local_min_right, local_max)
    maxmin = maximum([local_min_left, local_min_right])
    amplitude = abs( local_max - maxmin )
    return amplitude
end

function select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
    array_spikes_max_correct = Float64[]
    array_t_spikes_max_correct = Float64[]
    amplitudes = Float64[]
    for index in range(1, length(array_spikes_max)-1)
        amplitude = get_amplitude_spike(array_spikes_thresholds[index], array_spikes_thresholds[index+1], array_spikes_max[index])
        if amplitude >= 0.05
            push!(array_spikes_max_correct, array_spikes_max[index])
            push!(array_t_spikes_max_correct, array_t_spikes_max[index])
            push!(amplitudes, amplitude)
        end
    end
    return array_t_spikes_max_correct, array_spikes_max_correct, amplitudes
end
t_truncate(t) = floor(Int64, t / 2)