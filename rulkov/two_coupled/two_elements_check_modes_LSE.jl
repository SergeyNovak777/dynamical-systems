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

using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")


function two_coupled_rulkov_first_iteration(u, p)

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

    x1, y1, z1, x2, y2, z2 = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 = g2 * ( x_rp - x1 ) * xi(x2)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21, z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21)
    z1n = x1

    I12 = g1 * ( x_rp - x2 ) * xi(x1)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12, z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12)
    z2n = x2

    return SVector(x1n, y1n, z1n, I21, x2n, y2n, z2n, I12)
end

function two_coupled_rulkov(u, p, t)

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

    x1, y1, z1, I21prev, x2, y2, z2, I12prev = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 = γ_2 * I21prev + g2 * ( x_rp - x1 ) * xi(x2)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21, z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21)
    z1n = x1

    I12 = γ_1 * I12prev + g1 * ( x_rp - x2 ) * xi(x1)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12, z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12)
    z2n = x2

    return SVector{8}(x1n, y1n, z1n, I21, x2n, y2n, z2n, I12)
end

function get_params_two_coupled_rulkov()
    α = 4.6; σ =-0.1; μ = 0.001;
    β_syn = 0.0001; σ_syn = 1.0; k = 1.0;
    x_rp = -1.5; x_th = -0.8;
    γ_1 = 0.0; γ_2 = 0.0;
    g1 = 0.8; g2 = 0.2;
    return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2]
end

function pdf(array_peaks, threshold_range, ϵ)
    total_number_peaks = length(array_peaks)
    array_PDF = Float64[]
    for index in range(1, length(threshold_range), step = 1)
        number_peaks_crossing_theshold = length(findall( (threshold_range[index] + ϵ) .>= array_peaks .>= threshold_range[index] ) )
        PDF_for_theshold = number_peaks_crossing_theshold / total_number_peaks
        push!(array_PDF, PDF_for_theshold)
    end
    return array_PDF
end

function get_PDF(X, shift)
    total_count_X = length(X)
    array_PDF = Float64[]
    for index in 1:total_count_X
        count_X_i = count(X[index] -shift .<= X .<= X[index]+shift)
        PDF_X_i = count_X_i / total_count_X
        push!(array_PDF, PDF_X_i)
    end
    return array_PDF
end

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)
params = get_params_two_coupled_rulkov()
params[1] = 4.6;
params[2] = -0.1;
params = get_params_two_coupled_rulkov();
params[10] = 0.8;
params[11] = 0.2;
tspan = (0, 1000000);

u0 = [-0.1, -0.1, -0.1, -0.1, -0.1, -0.1];
u0_first_iteration = two_coupled_rulkov_first_iteration(u0, params);
prob = DiscreteProblem(two_coupled_rulkov, SVector{8}(u0_first_iteration), tspan, params);
sol = solve(prob);
t_start = 20000;
x_sum = sol[1, t_start:end] + sol[5, t_start:end]
t_range = sol.t[t_start:end]

ds = DeterministicIteratedMap(two_coupled_rulkov, sol[end], params)
Λs = lyapunov(ds, 500000)

t_plot_end = length(x_sum)
f = Figure(size = (1000, 400))
ax = Axis(f[1, 1])
lines!(ax, t_range[1:t_plot_end], x_sum[1:t_plot_end], linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)

#=
f = Figure(size = (1000, 400))
ax = Axis(f[1, 1])
lines!(ax, sol[1, t_start:t_plot_end], sol[5, t_start:t_plot_end], linewidth = 1.0, color = :blue)
lines!(ax, sol_asyn[1, t_start:t_plot_end], sol_asyn[5, t_start:t_plot_end], linewidth = 1.0, color = :red)
display(GLMakie.Screen(), f) =#

f = Figure(size = (400, 400))
ax = Axis3(f[1, 1])
scatter!(ax, sol[1, t_start:t_plot_end], sol[5, t_start:t_plot_end], sol[2, t_start:t_plot_end], markersize = 1.0, color = :blue)
display(GLMakie.Screen(), f)

println("syn: $Λs");