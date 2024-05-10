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
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/preprocessing_spike.jl")
end

using StaticArrays, DifferentialEquations, DynamicalSystems, Statistics, CairoMakie, GLMakie, JLD2

function clear_workspace()
    sol=nothing
    data_x1 = nothing
    data_x2 = nothing
    array_spikes_max = nothing
    array_t_spikes_max = nothing
    array_spikes_thresholds = nothing
    array_t_spikes_thresholds = nothing
    EE = nothing
    t_EE = nothing
    GC.gc()
end

parameters = FHN2_try3_params()
tspan = (0.0, 300000.0)
parameters[7] = 0.09
parameters[8] = 75.74
u0 = [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525] 
u0 = SVector{5}(u0)
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
sol = solve(prob, DP8(), adaptive = true, abstol = 1e-11, reltol = 1e-11, maxiters = 1e8);
len_sol = length(sol.t)
tstart = t_truncate(len_sol); tend = len_sol

data_x1 = [sol[1, :], sol.t]
sol=nothing; GC.gc(); sleep(3)

array_spikes_max, array_t_spikes_max = get_peaks(data_x1; level_zero = "nothing")
array_spikes_thresholds, array_t_spikes_thresholds = get_peaks_neg(data_x1)
baseline = Statistics.mean(array_spikes_thresholds)
println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")

if check_timeseries(array_spikes_max, array_spikes_thresholds) == true
    array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
else
    drop_false_start_end(array_t_spikes_max, array_spikes_max, array_t_spikes_thresholds)
    array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
end

println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")
println("count correct spike: $(length(array_spikes_max_correct))")

Hs_x  = Hs_above(amplitudes, 8)
index_EE = findall(x-> x >= Hs_x, array_spikes_max_correct)
t_EE = array_t_spikes_max_correct[index_EE]
EE = array_spikes_max_correct[index_EE]
len_EE = length(EE)
println("count EE: $(len_EE)")

IEI = Float64[]

for index in range(1, len_EE-1)
    t_EE_i = t_EE[index]
    t_EE_i_plus_one = t_EE[index+1]
    IEI_i = t_EE_i_plus_one - t_EE_i
    push!(IEI, IEI_i)
end
len_IEI = length(IEI)
IEI

#=t_plot_start = 1
t_shift = 150000
width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
lines!(ax, data_x1[2][t_plot_start:t_plot_start+t_shift], data_x1[1][t_plot_start:t_plot_start+t_shift], color = :blue, linewidth = 0.5)
scatter!(ax, t_EE, EE, markersize = 10, color = :deeppink)
xlims!(ax, data_x1[2][t_plot_start], data_x1[2][t_plot_start+t_shift])
hlines!(Hs_x, color = "red", linewidth = 2.0, linestyle = :dash)
hlines!(baseline, color = :blue, linewidth = 2.0, linestyle = :dash)
#IEI
lines!(ax, [t_EE[1], t_EE[2]], [0.0, 0.0], color = :black, linewidth = 1.0)
text!(ax, [t_EE[1], t_EE[2]], [0.0, 0.0], text = "IEI: $(IEI[1])")
lines!(ax, [t_EE[2], t_EE[3]], [0.0, 0.0], color = :black, linewidth = 1.0)
text!(ax, mean([t_EE[2], t_EE[3]]), 0.0, text = "IEI: $(IEI[2])")
display(GLMakie.Screen(), f)=#

width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"count", ylabel =  L"IEI")
lines!(ax, IEI, color = :black, linewidth = 1.0)
display(f)
