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
using StaticArrays, DynamicalSystems, DifferentialEquations, CairoMakie
include("/home/sergey/work/repo/dynamical-systems/system.jl")

function plot_phase_space(sol, tstart, tend, indexs, labels, path_to_save)
    f = Figure(size = (900, 600))
    ax = Axis3(f[1, 1], xlabel = labels[1], ylabel = labels[2], zlabel = labels[3], azimuth = -0.3pi,
        xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize, xticklabelsize= ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize)
    lines!(ax, sol[indexs[1], tstart:tend], sol[indexs[2], tstart:tend], sol[indexs[3], tstart:tend])
    display(f)
    save(path_to_save, f)
end

function plot_timeseries(X, trange)
    f = Figure(size = (1000, 250))
    ax = Axis(f[1,1], xlabel = L"time", ylabel = L"x_i")
    lines!(ax, trange, X[1], color = :red, linewidth = 1.5)
    lines!(ax, trange, X[2], color = :green, linewidth = 1.5)
    display(f)
end

function FHN_reduced(u, p, t)
    x1, y1, x2, y2 = u
    ϵ, a, g, k, σ, α, k1, k2 = p

    ϕ1 =  atan(y1, x1)
    ϕ2 =  atan(y2, x2)
    Iϕ2 = g / ( 1.0 + exp( k * ( cos(σ/2) - cos( ϕ2 - α - σ\2 ) ) ) )
    Iϕ1 = g / ( 1.0 + exp( k * ( cos(σ/2) - cos( ϕ1 - α - σ\2 ) ) ) )

    x1dt = x1 - x1^3 / 3 - y1 + Iϕ2 + ( k1 + k2 * ( y1 - y2 )^2 ) * (x2 - x1)
    x1dt = x1dt / ϵ
    y1dt = x1 - a 
    x2dt = x2 - x2^3 / 3 - y2 + Iϕ1 + ( k1 + k2 * (  y1 - y2 )^2 ) * (x1 - x2)
    x2dt = x2dt / ϵ
    y2dt = x2 - a

    return SVector{4}(x1dt, y1dt, x2dt, y2dt)
end

function get_params_FHN_reduced()
    ϵ = 0.01; a = -1.01;
    g = 0.1; k = 50.0; σ = 50.0 * pi / 180; α = 160.0 * pi / 180;
    k1 = 0.0; k2 = 0.0
    return [ ϵ, a, g, k, σ, α, k1, k2]
end

Ttr = 5000
t_LSE = 35000
integ_set = (alg = DP8(), adaptive = true, abstol = 1e-11, reltol = 1e-11)
parameters_fhn_reduced = get_params_FHN_reduced()
parameters_fhn_reduced[7] = 0.09
parameters_fhn_reduced[8] = 0.0
u0_fhn_reduced = [-2.0, 0.5, -1.0, 0.325]
t_solve = 1000
ds = CoupledODEs(FHN_reduced, u0_fhn_reduced, parameters_fhn_reduced, diffeq = integ_set)
traj, trange = trajectory(ds, t_solve)
X = [traj[:, 1], traj[:, 3]]
plot_timeseries(X, trange)