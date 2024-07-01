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

using StaticArrays, DifferentialEquations, JLD2, Statistics, CairoMakie, GLMakie

t_truncate(t) = floor(Int64, t / 2)
Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

alg = Vern9()
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
println("alg: $alg"); println("abstol: $abs_tol; reltol: $(rel_tol)")
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol,  max_iters = max_iters)

path_to_save = "/home/sergey/timeseries_k2_75_74_save_x1_x2/"
parameters = FHN2_try3_params()
parameters[7] = 0.09
parameters[8] = 75.74

u0_start = [-1.0537832727558796, -0.6375955063962165, -0.9461308916381171, -0.6265187262257061, -0.01107678017056557]; 
u0_start = SVector{5}(u0_start)

t_point = 100_000.0
tspan = (0.0, t_point)


prob = ODEProblem(FHN2_try3, u0_start, tspan, parameters)

sol = solve(prob, integrator_setting.alg, adaptive = true,
    abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
    maxiters = integrator_setting.max_iters);

len_sol = length(sol.t)
ttr = t_truncate(len_sol); tend = len_sol


data_x1 = [sol[1, ttr:tend], sol.t[ttr:tend]];
data_x2 = [sol[3, ttr:tend], sol.t[ttr:tend]];
 
data_local_max_x1 = get_local_max(data_x1)
data_local_min_x1 = get_local_min(data_x1)

data_local_max_x2 = get_local_max(data_x2)
data_local_min_x2 = get_local_min(data_x2)

drop_artifacts(data_local_max_x1, data_local_min_x1)
drop_artifacts(data_local_max_x2, data_local_min_x2)

all_amplitudes_x1 = get_amplitudes_all_events(data_local_max_x1[1], data_local_min_x1[1])
mean_amplitudes_x1 = Statistics.mean(all_amplitudes_x1)
println("mean amplitude x1: $mean_amplitudes_x1")
peaks_spikes_x1, t_peaks_spikes_x1, amplitudes_above_mean_x1 = select_spikes(data_local_min_x1[1], data_local_max_x1, mean_amplitudes_x1)

all_amplitudes_x2 = get_amplitudes_all_events(data_local_max_x2[1], data_local_min_x2[1])
mean_amplitudes_x2 = Statistics.mean(all_amplitudes_x2)
println("mean amplitude x2: $mean_amplitudes_x2")
peaks_spikes_x2, t_peaks_spikes_x2, amplitudes_above_mean_x2 = select_spikes(data_local_min_x2[1], data_local_max_x2, mean_amplitudes_x2)


Hs_x1 = Hs(amplitudes_above_mean_x1, 8);
Hs_x2 = Hs(amplitudes_above_mean_x2, 8);
println("Hs_x1: $Hs_x1");
println("Hs_x2: $Hs_x2");

index_EEs_x1 = findall(x-> x >= Hs_x1, peaks_spikes_x1)
peaks_EEs_x1 = peaks_spikes_x1[index_EEs_x1]
t_EEs_x1 = t_peaks_spikes_x1[index_EEs_x1]

index_EEs_x2 = findall(x-> x >= Hs_x2, peaks_spikes_x2)
peaks_EEs_x2 = peaks_spikes_x2[index_EEs_x2]
t_EEs_x2 = t_peaks_spikes_x2[index_EEs_x2]

println("len EE_x1: $(length(peaks_EEs_x1))");
println("len EE_x2: $(length(peaks_EEs_x2))");

array_IEI_x1 = get_IEI(t_EEs_x1)
array_PDF_IEI_x1 = get_PDF_IEI(array_IEI_x1; shift = 10)
Hs_IEI_coeff_8_x1 = Hs(array_IEI_x1, 8)
Hs_IEI_coeff_6_x1 = Hs(array_IEI_x1, 6)

array_IEI_x2 = get_IEI(t_EEs_x2)
array_PDF_IEI_x2 = get_PDF_IEI(array_IEI_x2; shift = 10)
Hs_IEI_coeff_8_x2 = Hs(array_IEI_x2, 8)
Hs_IEI_coeff_6_x2 = Hs(array_IEI_x2, 6)

labelsize = 40;
ticksize = 30;



#= 
f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI_{x1}", ylabel = L"PDF_{IEI_{x1}}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_IEI_x1, weights = array_PDF_IEI_x1, bins = 100)
vlines!(ax, Hs_IEI_coeff_8_x1, linewidth = 3.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_IEI_coeff_6_x2, linewidth = 3.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI_{x2}", ylabel = L"PDF_{IEI_{x2}}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_IEI_x2, weights = array_PDF_IEI_x2, bins = 100)
vlines!(ax, Hs_IEI_coeff_8_x2, linewidth = 3.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_IEI_coeff_6_x2, linewidth = 3.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f)


f = Figure()
ax = Axis(f[1, 1],
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
scatter!(ax, peaks_EEs_x1[2:end], peaks_EEs_x2[2:end])
display(GLMakie.Screen(), f) =#