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

using StaticArrays, DifferentialEquations, CairoMakie, GLMakie

function plot_timeseries(sol, tstart, tend, indexs, labels)
    indexx, indexy = indexs
    labelx, labely = labels
    f = Figure(size = (1000, 250))
    ax = Axis(f[1,1], xlabel = labelx, ylabel = labely)
    lines!(ax, sol.t[tstart:tend], sol[indexx, tstart:tend], color = :red)
    lines!(ax, sol.t[tstart:tend], sol[indexy, tstart:tend], color = :green)
    display(GLMakie.Screen(), f)
end

function plot_phase_space(sol, tstart, tend, indexs, labels)
    indexx, indexy, indexz = indexs
    labelx, labely, labelz = labels
    f = Figure(size = (400, 400))
    ax = Axis3(f[1,1], xlabel = labelx, ylabel = labely, zlabel = labelz)
    lines!(ax, sol[indexx, tstart:tend], sol[indexy, tstart:tend], sol[indexz, tstart:tend])
    display(GLMakie.Screen(), f)
end

function FHN2_try3(u, p ,t)
    x1, y1, x2, y2, z= u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))
    ρz = k1 + k2 * z ^ 2

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + ρz * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + ρz * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    dzdt = x1 - x2
    return SVector(dx1dt, dy1dt, dx2dt, dy2dt, dzdt)
end

function FHN2_try3_params()
    ϵ = 0.01; a = -1.01;
    g = 0.1; k = 50.0; σ = 50.0 * pi / 180; α = 160.0 * pi / 180;
    k1 = 0.0; k2 = 0.0
    return [ ϵ, a, g, k, σ, α, k1, k2]
end

parameters = FHN2_try3_params()
parameters[7] = 0.0
parameters[8] = 0.0
tspan = (0.0, 1000)

#u0 = SVector{5}([1.0, 0.0, 0.01, -1.0, 0.0])
u0 = SVector{5}([0.0, 1.0, 1.7, -1.0,0.0])
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
sol = solve(prob, RK4(), adaptive = false, dt = 0.001, maxiters = 5e6);

tstart = Int(length(sol) / 2); tend = length(sol)
indexx, indexy, indexz = 1, 3, 4
indexs = [indexx, indexy, indexz]
labelx, labely, labelz = "x1", "x2", "y2"
labels = [labelx, labely, labelz]


plot_phase_space(sol, tstart, tend, indexs, labels)
#plot_timeseries(sol, tstart, tend, [1,3], ["time", L"x_i"])