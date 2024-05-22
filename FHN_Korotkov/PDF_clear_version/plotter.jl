function plot_PDF_hist(array_IEI, array_PDF_IEI, bins_)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
    hist!(ax, array_IEI, weights = array_PDF_IEI, bins = bins_) #length(array_IEI))
    display(GLMakie.Screen(), f)
end

function plot_amplitudes(all_amplitudes, Hs_x; width_window = 1000, height_window = 300)
    f = Figure(size = (width_window, height_window))
    ax = Axis(f[1, 1], ylabel = L"amplitudes")
    lines!(ax, all_amplitudes)
    hlines!(ax, Hs_x, linewidth = 1.0, linestyle = :dash, color = :red)
    display(GLMakie.Screen(), f)
end

function plot_PDF_scatter(array_IEI, array_PDF_IEI; width_window = 1000, height_window = 300)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = L"IEI", ylabel = L"PDF_{IEI}", yscale = log10)
    scatter!(ax, array_IEI, array_PDF_IEI)
    display(GLMakie.Screen(), f)
end

function plot_timeseries(data, data_local_max, data_local_min,
    t_peaks_spikes, peaks_spikes, Hs_x,
    t_EEs, peaks_EEs,
    t_start, t_end;
    width_window = 1000, height_window = 300, clear_plot = false)
    f = Figure(size = (width_window, height_window))
    ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_1")
    lines!(ax, data[2][t_start:t_end], data[1][t_start:t_end])
    if clear_plot == false
        scatter!(ax, data_local_max[2], data_local_max[1], color = :green, markersize = 10)   
        scatter!(ax, t_peaks_spikes, peaks_spikes, color = :orange, markersize = 10)                                                                                                                                        
        scatter!(ax, data_local_min[2], data_local_min[1], color = :blue, markersize = 10)
        hlines!(ax, Hs_x, linewidth = 1.0, linestyle = :dash, color = :red)
        scatter!(ax, t_EEs, peaks_EEs, color = :deeppink, markersize = 10)
        xlims!(ax, data[2][t_start], data[2][t_end])
    end
    display(GLMakie.Screen(), f)
end