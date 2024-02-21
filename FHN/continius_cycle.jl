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
    include("/home/sergey/work/repo/dynamical-systems/system.jl")
end

using DifferentialEquations, DynamicalSystems, StaticArrays, CairoMakie

function plot_several_timeseries(sols, tstart, tend)
    fig = Figure(size=(800, 400))
    for i in 1:6
        ax = Axis(fig[i, 1], xticksvisible=false, yticksvisible=false,
                  xticklabelsvisible=false, yticklabelsvisible=false)
        lines!(ax, sols[i].t[tstart:tend], sols[i][1, tstart:tend], color = :red)
        lines!(ax, sols[i].t[tstart:tend], sols[i][3, tstart:tend], color = :green)
    end
    rowgap!(fig.layout, 0)
    display(fig)
end

function continuation_limit_cycles(sys, u0s, tspan, params, integ_setting)

    function prob_func(prob, index, repeat)
        remake(prob, u0 = u0s[index][end])
    end

    prob = ODEProblem(sys, u0s[1][end], tspan, params)
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(ensemble_prob, alg = integ_setting.alg, adaptive = integ_setting.adaptive, dt = integ_setting.dt, maxiters = 5e6,  trajectories = 6)
end

function continuation_limit_cycles_first(sys, u0s, tspan, params, integ_setting)

    function prob_func(prob, index, repeat)
        remake(prob, u0 = u0s[index, :])
    end

    prob = ODEProblem(sys, u0s[1, :], tspan, params)
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
    sim = solve(ensemble_prob, alg = integ_setting.alg, adaptive = integ_setting.adaptive, dt = integ_setting.dt, maxiters = 5e6, trajectories = 6)
end

integrator_setting = (alg = RK4(), adaptive = false, dt = 0.001)
params = two_coupled_fhn_get_params()
params[7] = 0.09
sys = two_coupled_fhn
tspan = (0.0, 1000.0)

#-------------------------------------------
u0s = zeros(6, 5)
u0s[1, :] = SVector(-1.0690237112785876, -0.5773325275931365, -1.1977568518748323, -0.5654141361024896, 1.60676442231909478)
u0s[2, :] = SVector(-1.1977568518748323, -0.5654141361024896, -1.0690237112785876, -0.5773325275931365, 2.50676442231909478)
u0s[3, :] = SVector(-1.1977568518748323, -0.5654141361024896, -1.0690237112785876, -0.5773325275931365, 2.00676442231909478)
u0s[4, :] = SVector(-1.1977568518748323, -0.5654141361024896, -1.0690237112785876, -0.5773325275931365, 1.00676442231909478)
u0s[5, :] = SVector(-1.1977568518748323, -0.5654141361024896, -1.0690237112785876, -0.5773325275931365, 0.50676442231909478)
u0s[6, :] = SVector(-1.1977568518748323, -0.5654141361024896, -1.0690237112785876, -0.5773325275931365, 0.00676442231909478)
u0s_static = SMatrix{6, 5}(u0s)
#-------------------------------------------

lengthrange = 50
k2range = range(0.0, 0.1, length = lengthrange)

tstart = 1; tend = 100000;

for index in range(1, lengthrange)
    params[8] = k2range[index]
    if index == 1
        global u0s_cont = continuation_limit_cycles_first(sys, u0s, tspan, params, integrator_setting)
        println("index: $(index)");flush(stdout)
        for j_index in range(1,6)
            println("last point from cycle-$(j_index): $(u0s_cont[j_index][end])")
        end
        println("---------------------");flush(stdout)
        println("");flush(stdout)
        plot_several_timeseries(u0s_cont, tstart, tend)
        continue
    end
    
    sim = continuation_limit_cycles(sys, u0s_cont, tspan, params, integrator_setting)
    u0s_cont = sim
    plot_several_timeseries(sim, tstart, tend)
    
    println("index: $(index)");flush(stdout)
    for j_index in range(1,6)
        println("last point from cycle-$(j_index): $(sim[j_index][1])")
    end
    println(">>>");flush(stdout)
    for j_index in range(1,6)
            println("last point from cycle-$(j_index): $(sim[j_index][end])")
    end
    println("---------------------");flush(stdout)
    println("");flush(stdout)
end
