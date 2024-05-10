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

path_to_save = "/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/sol_k2=75.74/"
iteration = 5
namefile_sol_x1 = "$(iteration)_sol_x1.jld"
namefile_sol_x2 = "$(iteration)_sol_x2.jld"
namefile_EE = "$(iteration)_EE.jld"
full_path_to_save_sol_x1 = path_to_save*namefile_sol_x1
full_path_to_save_sol_x2 = path_to_save*namefile_sol_x2
full_path_to_save_EE = path_to_save*namefile_EE

parameters = FHN2_try3_params()
t_point = 250000.0
tspan = (t_point * 4, t_point * 5)
parameters[7] = 0.09
parameters[8] = 75.74
u0 = [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525]
# 1 [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35]
# 2 [-1.0598231620372562, -0.6028280480069199, -0.9861681398046932, -0.6347932405747905, 0.03196519256785591]
# 3 [-1.0692095385745832, -0.6305096219263583, -0.9764407236291518, -0.6215810858043794, -0.008928536121969451]
# 4 [-0.9661682589719547, -0.6105027913464504, -1.027143212047144, -0.6372173689696131, 0.026714577623503304]
# 5 [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525]
# 6 
u0 = SVector{5}(u0)
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
sol = solve(prob, DP8(), adaptive = true, abstol = 1e-11, reltol = 1e-11, maxiters = 1e8);
len_sol = length(sol.t)
tstart = t_truncate(len_sol); tend = len_sol

data_x1 = [sol[1, :], sol.t]
data_x2 = [sol[3, :], sol.t]

array_spikes_max, array_t_spikes_max = get_peaks(data_x1; level_zero = "nothing")
array_spikes_thresholds, array_t_spikes_thresholds = get_peaks_neg(data_x1)
baseline = Statistics.mean(array_spikes_thresholds)
println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")

if check_timeseries(array_spikes_max, array_spikes_thresholds) == true
    # drop spike with small amplitude and save amplitude spike
    array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
else
    drop_false_start_end(array_t_spikes_max, array_t_spikes_thresholds)
    array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
end

println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")
println("count correct spike: $(length(array_spikes_max_correct))")

Hs_x  = Hs_above(amplitudes, 8)
index_EE = findall(x-> x >= Hs_x, array_spikes_max_correct)
t_EE = array_t_spikes_max_correct[index_EE]
EE = array_spikes_max_correct[index_EE]
println("count EE: $(length(EE))")

jldsave(full_path_to_save_sol_x1; data_x1)
jldsave(full_path_to_save_sol_x2; data_x2)
jldsave(full_path_to_save_EE, t_EE = t_EE, peaks_EE = EE)

#=width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
lines!(ax, sol.t[tstart:tstart+100000], sol[1, tstart:tstart+100000], color = :blue, linewidth = 0.5)
#scatter!(ax, t_EE, EE, markersize = 10, color = :deeppink)
scatter!(ax, array_t_spikes_max_correct, array_spikes_max_correct, markersize = 10, color = :green)
scatter!(ax, array_t_spikes_thresholds, array_spikes_thresholds, markersize = 10, color = :blue)
xlims!(ax, sol.t[tstart], sol.t[tstart+100000])
hlines!(Hs_x, color = "red", linewidth = 2.0, linestyle = :dash)
hlines!(baseline, color = :blue, linewidth = 2.0, linestyle = :dash)
display(GLMakie.Screen(), f)=#