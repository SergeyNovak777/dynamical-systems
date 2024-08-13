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
    include("/home/sergey/work/repo/dynamical-systems/system.jl");
end
using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie
using BenchmarkTools

parameters = get_params_three_coupled_rulkov()
parameters[10] = 1.0;
parameters[11] = 7.0;
tspan = (0, 300_000);

u0 = [-0.8653188666707976, -2.9297157423087055, -0.8883590191587176, -0.0, -0.0, -0.8462350817031328, -2.935527209875546, -0.8658876364860313, -0.0, -0.0, -0.9763100224806378, -2.9412769453476204, -0.9837730458388576, -0.0, -0.0]
u0_first_iteration = three_coupled_rulkov_first_iteration(u0, parameters);

prob_atr1 = DiscreteProblem(three_coupled_rulkov, u0_first_iteration, tspan, parameters);
sol_atr1 = solve(prob_atr1);


Ttr = 200_000;
point_from_attractor_atr1 = sol_atr1[:, end]

#= x_sum_atr1 = sol_atr1[1, Ttr:end] + sol_atr1[6, Ttr:end] + sol_atr1[11, Ttr:end]
sol_t_atr1 = sol_atr1.t[Ttr:end]
 =#
ds_atr1 = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor_atr1, parameters)

#Λs_atr1 = lyapunov(ds_atr1, 50_000)

#Λss_atr1 = lyapunovspectrum(ds_atr1, 50_000)

#@benchmark solve(prob_atr1);

@benchmark lyapunovspectrum(ds_atr1, 50_000)