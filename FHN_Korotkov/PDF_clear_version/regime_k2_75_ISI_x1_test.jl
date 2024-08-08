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
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl")
end

using StaticArrays, DifferentialEquations, JLD2, Statistics
using CairoMakie, GLMakie 

t_truncate(t) = floor(Int64, t / 2)
Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

alg = Vern9()
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
#= println("alg: $alg"); println("abstol: $abs_tol; reltol: $(rel_tol)")
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol, max_iters = max_iters) =#

path_to_save = "/home/sergey/timeseries_k2_75_74_save_x1_x2/"
parameters = FHN2_try3_params()
parameters[7] = 0.09
parameters[8] = 75.0;
println("k2: $(parameters[8])")
u0 = [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35]   
u0 = SVector{5}(u0)

t_point = 100_000.0;
tspan = (0.0, t_point);

prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = true,
abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
maxiters = integrator_setting.max_iters);

len_sol = length(sol.t)
ttr = t_truncate(len_sol); tend = len_sol
data_x1 = [sol[1, ttr:tend], sol.t[ttr:tend]]


data_local_max = get_local_max(data_x1)
data_local_min = get_local_min(data_x1)

drop_artifacts(data_local_max, data_local_min)

println("length local maxs: $(length(data_local_max[1]))")
println("length local mins: $(length(data_local_min[1]))")


all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
mean_amplitudes = Statistics.mean(all_amplitudes)
println("mean amplitude: $mean_amplitudes")
peaks_spikes, t_peaks_spikes, amplitudes_above_mean = select_spikes(data_local_min[1], data_local_max, mean_amplitudes) # Statistics.std(all_amplitudes))

data_local_max = nothing; data_local_min = nothing; GC.gc();

println("count spikes: $(length(peaks_spikes))")
Hs_x = Hs(amplitudes_above_mean, 8)

index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
peaks_EEs = peaks_spikes[index_EEs]
t_EEs = t_peaks_spikes[index_EEs]
println("count EEs: $(length(peaks_EEs))")

labelsize = 35;
ticksize = 20;

f = Figure()
ax = Axis(f[1, 1], xlabel = L"t", ylabel = L"x_1",
xgridvisible = false, ygridvisible = false,
xlabelsize = labelsize, ylabelsize = labelsize, xticklabelsize = ticksize, yticklabelsize = ticksize)
lines!(ax, sol.t[1:250_000], sol[3, 1:250_000], linewidth = 1.0)
display(GLMakie.Screen(), f)


f = Figure()
ax = Axis(f[1, 1], xlabel = L"t", ylabel = L"peaks_{x1}",
xgridvisible = false, ygridvisible = false,
xlabelsize = labelsize, ylabelsize = labelsize, xticklabelsize = ticksize, yticklabelsize = ticksize)
lines!(ax, t_peaks_spikes[1:end], peaks_spikes[1:end], linewidth = 1.0)
display(GLMakie.Screen(), f)


array_ISI = get_IEI(t_peaks_spikes);
array_PDF_ISI = get_PDF_IEI(array_ISI; shift = 10)

Hs_ISI_coeff_8 = Hs(array_ISI, 8)
Hs_ISI_coeff_6 = Hs(array_ISI, 6)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"ISI", ylabel = L"PDF_{ISI}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_ISI, weights = array_PDF_ISI, bins = 50)
vlines!(ax, Hs_ISI_coeff_8, linewidth = 5.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_ISI_coeff_6, linewidth = 5.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"t_{EE}", ylabel = L"IEI",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_peaks_spikes[2:end], array_ISI, linewidth = 1.0)
hlines!(ax, Hs_ISI_coeff_8, linewidth = 5.0, linestyle = :dash, color = :red)
hlines!(ax, Hs_ISI_coeff_6, linewidth = 5.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f)