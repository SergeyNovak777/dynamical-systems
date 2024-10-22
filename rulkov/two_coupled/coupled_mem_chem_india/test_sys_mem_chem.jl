username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

function rulkov_two_coupled_mem(u, p, t)
    x1, y1, x2, y2, ϕ = u;
    α, μ, k0, k1 ,k2, g21, g12, β_syn, σ_syn, x_rp, x_th = p

    function ξ(x, x_th)
        if x > x_th
            return 1.0;
        else
            return 0.0;
        end
    end

    x1_new = α / (1 + x1^2) + y1 + k0 * tanh(ϕ) * (x2 - x1) + g21 * (x_rp - x1) * ξ(x2, x_th) * β_syn / 2;
    y1_new = y1 - μ * (x1 + 1 + g21 * (x_rp - x1) * ξ(x2, x_th) * σ_syn / 2);
    x2_new = α / (1 + x2^2) + y2 + k0 * tanh(ϕ) * (x1 - x2) + g12 * (x_rp - x2) * ξ(x1, x_th) * β_syn / 2;
    y2_new = y2 - μ * (x2 + 1 + g12 * (x_rp - x2) * ξ(x1, x_th) * σ_syn / 2);
    ϕ_new = k1 * ϕ + k2 * (x1 - x2);

    return SVector(x1_new, y1_new, x2_new, y2_new, ϕ_new);
end


function get_params_rulkov_two_coupled_mem()
    α = 4.8; μ = 0.1; k0 = 1.5; k1 = 1.0; k2 = 0.01;
    g21 = 0.01; g12 = 0.01; x_th = -0.8; x_rp = -1.5;
    β_syn = 0.0001; σ_syn = 1.0;
    return [α, μ, k0, k1 ,k2, g21, g12, x_th, x_rp, β_syn, σ_syn];
end


params = get_params_rulkov_two_coupled_mem()
tspan = (0, 1500_000);

#=
α = 4.8;
u0 = [1.0, -2.0, 4.0, -1.1, 0.1];
=#

u0 = [1.0, -2.0, 4.0, -1.1, 0.1];
prob = DiscreteProblem(rulkov_two_coupled_mem, SVector{5}(u0), tspan, params);
sol = solve(prob);

ds = DeterministicIteratedMap(rulkov_two_coupled_mem, sol[end], params)
Λs = lyapunovspectrum(ds, 500_000)

t_start = 500_000;
t_plot_end = t_start + 700_000;

GLMakie.activate!();
f = Figure(size = (1000, 400))
ax = Axis(f[1, 1])
lines!(ax, sol.t[t_start:t_plot_end], sol[1, t_start:t_plot_end], linewidth = 1.0, color = :black)
display(GLMakie.Screen(), f)

f = Figure(size = (400, 400))
ax = Axis(f[1, 1])
lines!(ax, sol[1, t_start:t_plot_end], sol[2, t_start:t_plot_end], linewidth = 1.0, color = :black)
display(GLMakie.Screen(), f)

println("LLE: $Λs");