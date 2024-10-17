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

    function right_part_x(x, y, z, α)
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

    function right_part_y(x, y, sum_inhibitory, μ, σ)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x, x_th)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    ρ(L, k1 ,k2) = k1 + k2 * L^2;

    x1, y1, z1, x2, y2, z2, L = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, _, _, g1, g2, k1 ,k2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 = g2 * ( x_rp - x1 ) * xi(x2, x_th)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21 + ρ(L, k1, k2) * (x2 - x1), z1, α)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21, μ, σ)
    z1n = x1

    I12 = g1 * ( x_rp - x2 ) * xi(x1, x_th)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12 + ρ(L, k1 ,k2) * (x1 - x2), z2, α)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12, μ, σ)
    z2n = x2

    L = x1 - x2;

    return SVector(x1n, y1n, z1n, I21, x2n, y2n, z2n, I12, L)
end

function two_coupled_rulkov(u, p, t)

    function right_part_x(x, y, z, α)
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

    function right_part_y(x, y, sum_inhibitory, μ, σ)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x, x_th)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    ρ(L, k1 ,k2) = k1 + k2 * L^2;

    x1, y1, z1, I21prev, x2, y2, z2, I12prev, L = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2, k1, k2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 = γ_2 * I21prev + g2 * ( x_rp - x1 ) * xi(x2, x_th)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21 + ρ(L, k1, k2) * (x2 - x1), z1, α)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21, μ, σ)
    z1n = x1

    I12 = γ_1 * I12prev + g1 * ( x_rp - x2 ) * xi(x1, x_th)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12 + ρ(L, k1, k2) * (x1 - x2), z2, α)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12, μ, σ)
    z2n = x2

    L = x1 - x2;
    
    return SVector{9}(x1n, y1n, z1n, I21, x2n, y2n, z2n, I12,L)
end

function get_params_two_coupled_rulkov()
    α = 3.9; σ =1.0; μ = 0.001;
    β_syn = 0.0001; σ_syn = 1.0; k = 1.0;
    x_rp = -1.5; x_th = -0.8;
    γ_1 = 0.0; γ_2 = 0.0;
    g1 = 0.0; g2 = 0.0;
    k1 = 0.0; k2 = 0.0;
    return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2, k1, k2]
end
params = get_params_two_coupled_rulkov()
params[10] = 4.0;
params[11] = 2.5;
params[12] = 0.3;
params[13] = 0.005;
tspan = (0, 1000000);

u0 = [1.2, -0.3, 0.7, 2.1, 0.7, -0.1, 1.2 + 2.1];
u0_first_iteration = two_coupled_rulkov_first_iteration(u0, params);
prob = DiscreteProblem(two_coupled_rulkov, SVector{9}(u0_first_iteration), tspan, params);
sol = solve(prob);


t_start = 200_000;
x_sum = sol[1, t_start:end] + sol[5, t_start:end]
t_range = sol.t[t_start:end]

ds = DeterministicIteratedMap(two_coupled_rulkov, sol[end], params)
Λs = lyapunov(ds, 500000)

CairoMakie.activate!();

t_plot_end = t_start + 50_000;
f = Figure(size = (1000, 400))
ax = Axis(f[1, 1])
lines!(ax, t_range[1:t_plot_end], x_sum[1:t_plot_end], linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)

f = Figure(size = (400, 400))
ax = Axis3(f[1, 1])
scatter!(ax, sol[1, t_start:t_plot_end], sol[5, t_start:t_plot_end], sol[2, t_start:t_plot_end], markersize = 1.0, color = :blue)
display(GLMakie.Screen(), f)

println("LLE: $Λs");