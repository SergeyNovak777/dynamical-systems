# load model
include("/home/nova/work/repo_ds/dynamical-systems/tetrapartite synapse/mode_TM.jl");

username = "nova";
pathtorepo = "/home/nova/work/repo_ds/dynamical-systems";
using Pkg;
Pkg.activate(pathtorepo * "/env/integrate/");

using StaticArrays, DifferentialEquations, DynamicalSystems ;
using CairoMakie, GLMakie;

function get_par()
    τ = 0.013;  τD = 0.15; τF = 1.0;  τy = 1.8;
    α = 1.58; U0 = 0.23; ΔU0 = 0.305; I0 = -1.6;
    J = 3.07;  β = 0.4375;
    xth = 0.75; yth = 0.4;
    return [τ, τD, τF, τy, J, α, U0, ΔU0, I0, β, xth, yth];
end

function TM4glial_help()
    return "E-1, x-2, u-3, y-4 \n τ-1, τD-2, τF-3, τy-4, J-5, α-6, U0-7, ΔU0-8, I0-9, β-10, xth-11, yth-12";
end

function plot_timeseries(t, x, tstart, tend;
    width = 1000, height = 300, ylb = "x", lbsize = 35, ticksize = 25, grid = false, inter  = false)

    f= Figure(resolution = (width, height))
    axis  = Axis(f[1, 1], xlabel = L"time", ylabel = ylb,
    xlabelsize = lbsize, ylabelsize = lbsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = grid, ygridvisible = grid)

    lines!(axis, t[tstart:tend], x[tstart:tend], linewidth = 1.0, color = :black)
    if inter == false
        display(f)
    else
        display(GLMakie.Screen(), f)
    end
end

function main()

    par = get_par()
    init_cond = [3.3403239669724387, 0.1, 0.1, 0.03677942307955071];
    t_step = 0.001;

    integ_set = (alg = RK4(), adaptive = false, dt = t_step);

    ds = CoupledODEs(TM4glial, init_cond, par, diffeq = integ_set);
    
    return ds;

end