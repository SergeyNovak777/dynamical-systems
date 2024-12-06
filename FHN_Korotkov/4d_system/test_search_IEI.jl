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

using StaticArrays, DifferentialEquations, BenchmarkTools, JLD2, Statistics, CairoMakie, GLMakie

t_truncate(t) = floor(Int64, t / 2);
Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x);

alg = Vern9();
abs_tol = 1e-7;
rel_tol = 1e-7;
max_iters = 1e8;
#println("alg: $alg"); println("abstol: $abs_tol; reltol: $(rel_tol)")
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol,  max_iters = max_iters);

parameters = FHN2_try3_params();
parameters[7] = 0.09;
parameters[8] = 75.74;

u0_start = [1.7, 0.7, -1.4, 0.35]; 
u0_start = SVector{4}(u0_start);

t_point = 3_000_000;
tspan = (0.0, t_point);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters);

sol = solve(prob, integrator_setting.alg, adaptive = true,
    abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
    maxiters = integrator_setting.max_iters, save_idxs = [1], dense = false);

len_sol = length(sol.t)
ttr = t_truncate(len_sol); tend = len_sol

data_x1 = [sol[1, ttr:tend], sol.t[ttr:tend]];
print("start GC");
sol = nothing; GC.gc(); GC.gc();
data_local_max_x1 = get_local_max(data_x1)
data_local_min_x1 = get_local_min(data_x1)

drop_artifacts(data_local_max_x1, data_local_min_x1)

all_amplitudes_x1 = get_amplitudes_all_events(data_local_max_x1[1], data_local_min_x1[1])
mean_amplitudes_x1 = Statistics.mean(all_amplitudes_x1)
println("mean amplitude x1: $mean_amplitudes_x1")
peaks_spikes_x1, t_peaks_spikes_x1, amplitudes_above_mean_x1 = select_spikes(data_local_min_x1[1], data_local_max_x1, mean_amplitudes_x1)

Hs_x1 = Hs(amplitudes_above_mean_x1, 8);

println("Hs_x1: $Hs_x1");

index_EEs_x1 = findall(x-> x >= Hs_x1, peaks_spikes_x1)
peaks_EEs_x1 = peaks_spikes_x1[index_EEs_x1]
t_EEs_x1 = t_peaks_spikes_x1[index_EEs_x1]

println("len EE_x1: $(length(peaks_EEs_x1))");

array_IEI_x1 = get_IEI(t_EEs_x1)
array_PDF_IEI_x1 = get_PDF_IEI(array_IEI_x1; shift = 10)
Hs_IEI_coeff_8_x1 = Hs(array_IEI_x1, 8)
Hs_IEI_coeff_6_x1 = Hs(array_IEI_x1, 6)
println("Hs_IEI_coeff_8: $Hs_IEI_coeff_8_x1");
println("Hs_IEI_coeff_6: $Hs_IEI_coeff_6_x1");
labelsize = 40;
ticksize = 30;

f = Figure()
ax = Axis(f[1, 1], xlabel = L"IEI_{x1}", ylabel = L"PDF_{IEI_{x1}}", yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
hist!(ax, array_IEI_x1, weights = array_PDF_IEI_x1, bins = 100)
vlines!(ax, Hs_IEI_coeff_8_x1, linewidth = 3.0, linestyle = :dash, color = :red)
vlines!(ax, Hs_IEI_coeff_6_x1, linewidth = 3.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"t_{EE}", ylabel = L"IEI",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)
lines!(ax, t_EEs_x1[2:end], array_IEI_x1, linewidth = 1.0)
hlines!(ax, Hs_IEI_coeff_8_x1, linewidth = 5.0, linestyle = :dash, color = :red)
hlines!(ax, Hs_IEI_coeff_6_x1, linewidth = 5.0, linestyle = :dash, color = :green)
display(GLMakie.Screen(), f)