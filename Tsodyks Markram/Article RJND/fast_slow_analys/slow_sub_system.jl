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

using StaticArrays, DifferentialEquations, DynamicalSystems
xlims, ylims = 0.0 .. 1.0, 0.0 .. 1.0
box = xlims × ylims 
using CairoMakie

function slow_sub_system(u, p, t)
    x, y = u
    τD, U0, ΔU0, ythr, τy, β, xthr, E = p
    Uy = U0 + ΔU0 / ( 1 + exp( -50 * ( y - ythr ) ) )
    σy = 1 / ( 1 + exp( -20 * ( x - xthr ) ) )
    dxdt = ( 1 - x ) / τD - Uy * x * E
    dydt = -y / τy + β * σy
    return SVector{2}(dxdt, dydt)
end

function jac_slow_sub_system(u, p, t)

    x, y = u
    τD, U0, ΔU0, ythr, τy, β, xthr, E = p
    Uy = U0 + ΔU0 / ( 1 + exp( -50 * ( y - ythr ) ) )
    σy = 1 / ( 1 + exp( -20 * ( x - xthr ) ) )

    exponent_with_y = exp( -50 * ( y - ythr ) )
    exponent_with_x = exp( -20 * ( x - xthr ) )

    xx = -1 / τD - Uy * E
    xy = ( -50 * x * E * ΔU0 * exponent_with_y ) / ( 1 + exponent_with_y ) ^2
    yx = 20 * β * exponent_with_x / ( 1 + exponent_with_x )^2
    yy = -1/τy
    SMatrix{2,2}(xx, yx,
                xy, yy)
end
    

function get_parameters_slow_sub_system()
    τD = 0.07993; U0 = 0.265; ΔU0 = 0.305; ythr = 0.4; τy = 3.3; β = 0.3; xthr = 0.75; E = 15.0;
    return τD, U0, ΔU0, ythr, τy, β, xthr, E
end

function solver(prob, integrator_setting)
    if integrator_setting.adaptive == true
        solution = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol)
    else
        solution = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        dt = integrator_setting.dt,
        maxiters = integrator_setting.maxiters)
    end
    return solution
end

function plot_phase_space(X, xlbl, ylbl; window_width = 500, window_height = 500)
    x1, x2 = X
    fig = Figure(size = (window_width, window_height))
    ax = Axis(fig[1, 1], xlabel = xlbl, ylabel = ylbl)
    lines!(ax, x1, x2)
    display(fig)
end

sys = slow_sub_system
parameters = get_parameters_slow_sub_system()
u0 = SVector{2}(0.15, 0.23)
integrator_setting = (alg = Vern9(), adaptive = true, abstol = 1e-9, reltol = 1e-9)
tspan = (0.0, 1000.0)

ode_ds = CoupledODEs(sys, u0, parameters, diffeq = integrator_setting)
sol, trange = trajectory(ode_ds, 1000)
X = [ sol[:,1], sol[:, 2] ]

plot_phase_space(X, "x", "y")

fp, eigens, _ = fixedpoints(ode_ds, box, jac_slow_sub_system)