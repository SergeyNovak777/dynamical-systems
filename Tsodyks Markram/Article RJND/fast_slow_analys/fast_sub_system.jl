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
E_box = 0..30

function fast_sub_sys_ds(u, p ,t)
    Uy = p[4] + p[5] / ( 1.0 + exp( -50.0 * ( p[9] - p[6] ) ) )
    dEdt = -u[1] + p[2] * log( 1.0 + exp( (p[3] * Uy * p[8] * u[1] + p[7]) / p[2] ) )
    dEdt = dEdt / p[1]
    return SVector{1}(dEdt)
end

function jac_fast_sub_sys(u, p, t)
    Uy = p[4] + p[5] / ( 1.0 + exp( -50.0 * ( p[9] - p[6] ) ) )
    expression_under_exp = (p[3] * Uy * p[8] * u[1] + p[7]) / p[8]
    the_exp = exp(expression_under_exp)
    EE = -1.0 + the_exp * ( p[3] * Uy * p[8] + p[7] ) / ( 1.0 + the_exp )
    EE = EE / p[1]
    return SMatrix{1,1}(EE)
end

function get_param_fast_sub_sys()
    τ = 0.013; α = 1.58; J = 3.07; U0 = 0.265; ΔU0 = 0.305; ythr = 0.4; I0 = -1.6;
    x = 0.0; y = 0.0;
    return [ τ, α, J, U0, ΔU0, ythr, I0, x, y]
end

function get_help_fast_sub_sys(params)
    indexparams = "τ - 1, α - 2, J - 3, U0 - 4,
    ΔU0 - 5, ythr - 6, I0 - 7, x - 8, y - 9";
    nameparams = "τ, α, J, U0, ΔU0, ythr, I0, x, y";
    keyp = split(nameparams, ", ");
    dict = Dict(zip(keyp, params));
    return dict, indexparams;
end

parameters = get_param_fast_sub_sys()
u0 = 3.0
tspan = (0.0, 1000.0)
integrator_setting = (alg = Vern9(), adaptive = true, abstol = 1e-9, reltol = 1e-9)

ode_ds = CoupledODEs(fast_sub_sys_ds, u0, parameters, diffeq = integrator_setting)

jac_fast_sub_sys(u0, parameters, 0)

fixedpoints(ode_ds, E_box, jac_fast_sub_sys);


#=
using CairoMakie
xs, trange = trajectory(ode_ds, 1000)
sys = fast_sub_sys
function fast_sub_sys(u, p ,t)
    Uy = p[4] + p[5] / ( 1.0 + exp( -50.0 * ( p[9] - p[6] ) ) )
    dEdt = -u + p[2] * log( 1.0 + exp( (p[3] * Uy * p[8] * u + p[7]) / p[2] ) )
    dEdt = dEdt / p[1]
    return dEdt
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

function plot_timeseries(tspan, x, ylbl; window_width = 1000, window_height = 300)
    fig = Figure(size = ( window_width, window_height ))
    ax = Axis(fig[1, 1], xlabel = L"time", ylabel = ylbl)
    lines!(ax, tspan, x)
    display(fig)
end=#

#ode_de = ODEProblem(sys, u0, tspan, parameters)
#sol = solver(ode_de, integrator_setting)


#plot_timeseries(sol.t, sol.u, L"E")
#plot_timeseries(trange, xs[:,1], L"E")