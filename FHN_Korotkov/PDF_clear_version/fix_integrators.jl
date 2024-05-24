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

using StaticArrays, DifferentialEquations, Plots

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
tspan = (0.0, 200_000.0)
parameters[7] = 0.09
parameters[8] = 75.74
u0 = [-0.9859005363852416, -0.635253572091177, -1.0345181027025165, -0.636382088705782, 0.0011285166148596525] 
u0 = SVector{5}(u0)
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)

alg = Vern9();
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
integrator_setting = (alg = alg, abs_tol = abs_tol, rel_tol = rel_tol, max_iters = max_iters, dense = false)

sol = solve(prob, integrator_setting.alg, adaptive = true, abstol = integrator_setting.abs_tol, reltol = integrator_setting.rel_tol,
maxiters = integrator_setting.max_iters, dense = false);

len_sol = length(sol);
tstart = len_sol - 500;
tend = len_sol;
plot(sol.t[tstart:tend], sol[1, tstart:tend], fmt = :pdf)
ylims!(-1.1, -0.90)


#= f = Figure()
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_1")
lines!(ax, sol.t[tstart:tend], sol[1, tstart:tend], linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)  =#

#= tstart = tspan[1]; tend = tspan[2]
trange = range(tstart, tend, length = 5000000)
values_interp = sol(trange)[1, :]

sol = nothing; GC.gc() =#