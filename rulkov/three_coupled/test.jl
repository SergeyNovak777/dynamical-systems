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
end
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl")
using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie, JLD2

function three_coupled_rulkov_first_iteration(u, p)

    function right_part_x(x, y, z)
        if x <= 0.0
            return α/(1.0 - x) + y
        elseif (0.0 < x < (α + y)) && (z <=0.0)
            return α  +y
        elseif (x >= α + y) || (z > 0)
            return -1.0
        else
            return -1.0;
        end
    end

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, x2, y2, z2, x3, y3, z3 = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 = g2 * ( x_rp - x1 ) * xi(x2)
    I31 = g1 * ( x_rp - x1 ) * xi(x3)

    I12 = g1 * ( x_rp - x2 ) * xi(x1)
    I32 = g2 * ( x_rp - x2 ) * xi(x3)

    I13 = g2 * ( x_rp - x3 ) * xi(x1)
    I23 = g1 * ( x_rp - x3 ) * xi(x2)

    x1n = right_part_x(x1, y1 + (β_syn/k) * (I21 + I31), z1)    
    y1n = right_part_y(x1, y1, (σ_syn/k) * (I21 + I31))
    z1n = x1

    x2n = right_part_x(x2, y2 + (β_syn/k) * (I12 + I32), z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * (I12 + I32))
    z2n = x2

    x3n =  right_part_x(x3, y3 + (β_syn/k) * (I13 + I23), z3)
    y3n = right_part_y(x3, y3, (σ_syn/k) * (I13 + I23))
    z3n = x3

    return SVector(x1n, y1n, z1n, I21, I31, x2n, y2n, z2n, I12, I32, x3n, y3n, z3n, I13, I23)
end

function three_coupled_rulkov(u, p, t)

    function right_part_x(x, y, z)
        if x <= 0.0
            return α/(1.0 - x) + y
        elseif (0.0 < x < (α + y)) && (z <=0.0)
            return α  +y
        elseif (x >= α + y) || (z > 0)
            return -1.0
        else
            return -1.0;
        end
    end

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, I21prev, I31prev, x2, y2, z2, I12prev, I32prev, x3, y3, z3, I13prev, I23prev = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 =  γ_2 * I21prev + g2 * ( x_rp - x1 ) * xi(x2)
    I31 =  γ_1 * I31prev + g1 * ( x_rp - x1 ) * xi(x3)
    x1n = right_part_x(x1, y1 + (β_syn/k) * (I21 + I31), z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * (I21 + I31))
    z1n = x1

    I12 = γ_1 * I12prev + g1 * ( x_rp - x2 ) * xi(x1)
    I32 = γ_2 * I32prev + g2 * ( x_rp - x2 ) * xi(x3)
    x2n = right_part_x(x2, y2 + (β_syn/k) * (I12 + I32), z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * (I12 + I32))
    z2n = x2

    I13 = γ_2 * I13prev +  g2 * ( x_rp - x3 ) * xi(x1)
    I23 = γ_1 * I23prev + g1 * ( x_rp - x3 ) * xi(x2)
    x3n =  right_part_x(x3, y3 + (β_syn/k) * (I13 + I23), z3)
    y3n = right_part_y(x3, y3, (σ_syn/k) * (I13 + I23))
    z3n = x3

    return SVector(x1n, y1n, z1n, I21, I31, x2n, y2n, z2n, I12, I32, x3n, y3n, z3n, I13, I23)
end

function get_params_three_coupled_rulkov()
    α = 3.9; σ =1.0; μ = 0.001;
    β_syn = 0.0001; σ_syn = 1.0;
    x_rp = -1.5; x_th = -0.8;
    γ_1 = 0.0; γ_2 = 0.0;
    g1 = 1.0; g2 = 7.0;
    return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2]
end

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function CALCPDF_debil(spikes, threshold, ϵ)
    ee_counter = [sum(i->s<=i<s+ϵ, spikes) for s in threshold]
    pdf = ee_counter ./ length(spikes)
    return pdf
end

params = get_params_three_coupled_rulkov()
params[10] = 4.7;
params[11] = 5.0;
tspan = (0, 5000000);

u0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
u0_first_iteration = three_coupled_rulkov_first_iteration(u0, params);
prob = DiscreteProblem(three_coupled_rulkov, SVector{length(u0_first_iteration)}(u0_first_iteration), tspan, params);
sol = solve(prob);
Ttr = 1_000_000;
point_from_attractor = sol[:, Ttr]
x_sum = sol[1, :] + sol[6, :] + sol[11, :]

#= ds = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor, params)
Λs = lyapunov(ds, 500000)
println("Λs: $Λs") =#

data = [x_sum[Ttr:end], sol.t[Ttr:end]]
sol = nothing;
x_sum = nothing;

data_local_max = get_local_max(data)
data_local_min = get_local_min(data)

drop_artifacts(data_local_max, data_local_min)

Hs_xsum = Hs(data_local_max[1] ,6);
thesholds = range(0.1, 3, 5_000_000);
PDF_old = CALCPDF_debil(data_local_max[1], thesholds, 0.05);

EE_mapcopy = PDF_old;
EE_mapcopy = [ iszero(x) ? NaN : x for x in EE_mapcopy ];

labelsize = 40;
ticksize = 30;
CairoMakie.activate!();

path_to_save_timeseries = "/home/sergey/MEGA/dynamical-systems/Rulkov/Images/timeseries/";
path_to_save_PDF = "/home/sergey/MEGA/dynamical-systems/Rulkov/Images/PDF/";
path_to_save_data = "/home/sergey/MEGA/dynamical-systems/Rulkov/data/PDF/";

# timeseries of xsum
t_plot_start = 1; t_plot_end = 200_000;
f = Figure(size = (1000, 400));
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_{sum}",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(ax, data[2][t_plot_start:t_plot_end], data[1][t_plot_start:t_plot_end], linewidth = 1.0, color = :blue);
hlines!(ax, Hs_xsum, linestyle = :dash, color = :red, linewidth = 3.0);
display(GLMakie.Screen(), f);
save(path_to_save_timeseries*"timeseries__xsum_g1=4.7_g2=5.0_.eps", f)

# pdf old version
f = Figure();
ax = Axis(f[1, 1], xlabel = L"x_{sum}", ylabel = L"PDF",
xscale = log10, yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(thesholds, EE_mapcopy, linewidth = 1.0, color = :blue);
vlines!(ax, Hs_xsum, linewidth = 3.0, linestyle = :dash, color = :red);
display(GLMakie.Screen(), f);
save(path_to_save_PDF*"PDF_g1=4.7_g2=5.0_.eps", f)


#= jldsave(path_to_save_data*"PDF_g1=4.7_g2=5.0_jld2"; EE_mapcopy);
jldsave(path_to_save_data*"data_g1=4.7_g2=5.0_jld2"; data); =#


amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1]);
PDF_amplitudes = get_PDF(amplitudes, 0.05);
amplitudes_sorted = sort(amplitudes);
PDF_amplitudes_sorted = get_PDF(amplitudes_sorted, 0.05);
Hs_xsum_ampl = Hs(amplitudes, 6);

array_local_max_above_thr, array_t_local_max_above_thr, amplitudes_above_thr = select_spikes(data_local_min[1], data_local_max, 0.5);
println("minumum amplitude: $(minimum(amplitudes_above_thr))");
amplitudes_above_thr_sorted = sort(amplitudes_above_thr);
PDF_amplitudes_thr = get_PDF(amplitudes_above_thr_sorted, 0.05);
Hs_xsum_ampl_above_thr = Hs(amplitudes_above_thr_sorted, 6);

# pdf by amplitudes
CairoMakie.activate!();
f = Figure();
ax = Axis(f[1, 1], xlabel = L"amplitudes", ylabel = L"PDF",
xscale = log10, yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(amplitudes_sorted, PDF_amplitudes_sorted, linewidth = 1.0, color = :blue);
vlines!(ax, Hs_xsum_ampl, linewidth = 3.0, linestyle = :dash, color = :red);
display(GLMakie.Screen(), f);
save(path_to_save_PDF*"PDF_ampl_g1=4.7_g2=5.0_.eps", f);

# pdf by amplitudes above thr
CairoMakie.activate!();
f = Figure();
ax = Axis(f[1, 1], xlabel = "amplitudes_above_thr", ylabel = L"PDF",
xscale = log10, yscale = log10,
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(amplitudes_above_thr_sorted, PDF_amplitudes_thr, linewidth = 1.0, color = :blue);
vlines!(ax, Hs_xsum_ampl_above_thr, linewidth = 3.0, linestyle = :dash, color = :red);
display(GLMakie.Screen(), f);
save(path_to_save_PDF*"PDF_ampl_above_thr_g1=4.7_g2=5.0_.eps", f);

# timeseries of amplitudes
t_plot_start = 1; t_plot_end = 200_000;
f = Figure(size = (1000, 400));
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"amplitudes",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false);
lines!(ax, amplitudes[t_plot_start:t_plot_end], linewidth = 1.0, color = :blue);
hlines!(ax, Hs_xsum_ampl, linestyle = :dash, color = :red, linewidth = 3.0);
display(GLMakie.Screen(), f);
save(path_to_save_timeseries*"timeseries_ampl_g1=4.7_g2=5.0_.eps", f)