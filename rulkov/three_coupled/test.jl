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

params = get_params_three_coupled_rulkov()
params[10] = 3.1;
params[11] = 3.0;
tspan = (0, 1000000);

u0 = [1.0, 0.5, 0.2, -1.0, -0.5, -0.2, 1.5, 1.0, 0.4]
u0_first_iteration = three_coupled_rulkov_first_iteration(u0, params);
prob = DiscreteProblem(three_coupled_rulkov, SVector{length(u0_first_iteration)}(u0_first_iteration), tspan, params);
sol = solve(prob);
x_sum = sol[1, :] + sol[6, :] + sol[11, :]

ds = DeterministicIteratedMap(three_coupled_rulkov, sol[end], params)
Λs = lyapunov(ds, 500000)
println("Λs: $Λs")

t_start_plot = 300_000
t_end_plot = 300_800

f = Figure(size = (400, 400))
ax1 = Axis(f[1, 1])
ax2 = Axis(f[2, 1])
ax3 = Axis(f[3, 1])
lines!(ax1, sol.t[t_start_plot:t_end_plot], sol[1, t_start_plot:t_end_plot], linewidth = 1.0, color = :green)
lines!(ax2, sol.t[t_start_plot:t_end_plot], sol[6, t_start_plot:t_end_plot], linewidth = 1.0, color = :red)
lines!(ax3, sol.t[t_start_plot:t_end_plot], sol[11, t_start_plot:t_end_plot], linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)


t_start_plot = 300_000
t_end_plot = 1000_000
f = Figure(size = (400, 400))
ax = Axis3(f[1, 1])
scatter!(ax, sol[1, t_start_plot:t_end_plot], sol[6, t_start_plot:t_end_plot], sol[11, t_start_plot:t_end_plot], markersize = 1.0, color = :blue)
display(GLMakie.Screen(), f)