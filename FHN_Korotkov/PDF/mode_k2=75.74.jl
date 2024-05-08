if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
    include("/home/sergey/work/repo/dynamical-systems/system.jl")
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/pdf_function.jl")
end


using StaticArrays, DifferentialEquations, DynamicalSystems, Statistics, CairoMakie, GLMakie

t_truncate(t) = floor(Int64, t / 2)

parameters = FHN2_try3_params()
tspan = (0.0, 300000.0)
parameters[7] = 0.09
parameters[8] = 75.74
u0 = [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35]
u0 = SVector{5}(u0)
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
sol = solve(prob, DP8(), adaptive = true, abstol = 1e-11, reltol = 1e-11, maxiters = 1e8);
len_sol = length(sol.t)
tstart = t_truncate(len_sol); tend = len_sol

data = [sol[1, tstart:tend], sol.t[tstart:tend]]
array_spikes_max, array_t_spikes_max = get_peaks(data; level_zero = "nothing")
array_spikes_thresholds, array_t_spikes_thresholds = get_peaks_neg(data)
baseline = Statistics.mean(array_spikes_thresholds)
println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")

function check_timeseries(array_spikes_max, array_spikes_thresholds)
    if length(array_spikes_max)+1 == length(array_spikes_thresholds)
        flag = true
    else
        flag = false
    end
    return flag
end

function drop_false_start_end(array_t_spikes_max, array_t_spikes_thresholds)
    if array_t_spikes_max[1] < array_t_spikes_thresholds[1]
        popfirst!(array_t_spikes_max)
        popfirst!(array_spikes_max)
    elseif array_t_spikes_max[end] > array_t_spikes_thresholds[end]
        pop!(array_t_spikes_max)
        pop!(array_spikes_max)
    end
end

if check_timeseries(array_spikes_max, array_spikes_thresholds) == true
    # drop spike with small amplitude and save amplitude spike
else
    drop_false_start_end(array_t_spikes_max, array_t_spikes_thresholds)
end

function get_amplitude_spike(local_min_left, local_min_right, local_max)
    maxmin = maximum([local_min_left, local_min_right])
    amplitude = abs( local_max - maxmin )
    return amplitude
end

array_spikes_max_correct = Float64[]
array_t_spikes_max_correct = Float64[]
amplitudes = Float64[]
for index in range(1, length(array_spikes_max)-1)
    amplitude = get_amplitude_spike(array_spikes_thresholds[index], array_spikes_thresholds[index+1], array_spikes_max[index])
    if amplitude >= 0.1
       push!(array_spikes_max_correct, array_spikes_max[index])
       push!(array_t_spikes_max_correct, array_t_spikes_max[index])
       push!(amplitudes, amplitude)
    end
end

println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")
println("count correct spike: $(length(array_spikes_max_correct))")

minspike = minimum(amplitudes)
maxspike = maximum(amplitudes)

Hs_x  = Hs_above(amplitudes, 8)

width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
lines!(ax, sol.t[tstart:tstart+100000], sol[1, tstart:tstart+100000], color = :blue, linewidth = 0.5)
scatter!(ax, array_t_spikes_max_correct, array_spikes_max_correct, markersize = 10, color = :green)
scatter!(ax, array_t_spikes_thresholds, array_spikes_thresholds, markersize = 10, color = :blue)
xlims!(ax, sol.t[tstart], sol.t[tstart+100000])
hlines!(Hs_x, color = "red", linewidth = 2.0, linestyle = :dash)
hlines!(baseline, color = :blue, linewidth = 2.0, linestyle = :dash)
display(GLMakie.Screen(), f)
,,,