if Sys.iswindows()
    cd("D:\\work\\dynamical-systems\\env\\")
    using Pkg
    Pkg.activate("integrate")
    include("D:\\work\\dynamical-systems\\system.jl");
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
    include(pathtorepo*"/system.jl")
    
end

using DifferentialEquations, DynamicalSystems, CSV, DataFrames, DelimitedFiles

function preproc_df()
    
    df = DataFrame(CSV.File("/home/sergey/work/repo/dynamical-systems/rate model/randomizer/50000_u0s.csv"))
    df_filter = filter(r -> r.LE1 >= 0, df)

    u0s = df_filter[:,[:"sE_start",:"sI_start", :"rE_start", :"rI_start", :"Y_start"]]

    u0s = Matrix{Float64}(u0s) 
    return u0s
end

function init_ODE_prob(sys, params, u0_lc, time_attract)
    tspan = (0.0, time_attract)
    prob = ODEProblem(sys, SVector{length(u0_lc)}(u0_lc), tspan, params)
    return prob
end

function init_Coupled_ODE(sys, params, u0_lc, integ_set)
    ds = CoupledODEs(sys, u0_lc, params, diffeq = integ_set)
    return ds
end

function goto_attractor(prob, integrator_setting)

    if integrator_setting.adaptive == true
        traj = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol,
        save_everystep = false, save_start = false)
    else
        traj = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        dt = integrator_setting.dt, maxiters = integrator_setting.maxiters,
        save_everystep = false, save_start = false)
    end
    return traj[end]
end

function calculate_LSE(ds, time_calculate_LSE)
    LSE = lyapunovspectrum(ds, time_calculate_LSE)
    return LSE
end

function main()

    data, header = readdlm("/home/sergey/work/repo/dynamical-systems/rate model/randomizer/output.csv", ',', header=true);
    info_traj = DataFrame(data, vec(header))

    system = rate_model
    parameters = rate_model_get_params()
    t_attract = 1000.0
    t_LSE = 20000
    integrator_setting = (alg = Vern9(), adaptive = true, abstol = 1e-11, reltol = 1e-11)
    
    u0s = preproc_df()

    for index in range(1, length(u0s), step = 1)
        
        prob = init_ODE_prob(system, parameters, u0s[index, :], t_attract)
        last_point = goto_attractor(prob, integrator_setting)
        ds = init_Coupled_ODE(system, parameters, last_point, integrator_setting)
        LSE = calculate_LSE(ds, t_LSE)

        sE_sv, sI_sv, rE_sv, rI_sv, Y_sv = u0s[index, :]
        sE_ev, sI_ev, rE_ev, rI_ev, Y_ev = last_point
        l1, l2, l3, l4, l5 = LSE
        
        if index == 1
            global df = DataFrame("sE_start" => sE_sv, "sI_start" => sI_sv, "rE_start" => rE_sv, "rI_start" => rI_sv, "Y_start" => Y_sv,
                "sE_end" => sE_ev, "sI_end" => sI_ev, "rE_end" => rE_ev, "rI_end" => rI_ev, "Y_end" => Y_ev,
                "LE1" => l1, "LE2" => l2, "LE3" => l3, "LE4" => l4, "LE5" => l5)
        else
            push!(df, (sE_sv, sI_sv, rE_sv, rI_sv, Y_sv, sE_ev, sI_ev, rE_ev, rI_ev, Y_ev, l1, l2, l3, l4, l5))
        end

        if mod(index, 10) == 0
            pathtofile = "/home/sergey/work/repo/dynamical-systems/rate model/randomizer/"
            namefile = "continue_u0s.csv"
            fullpath = pathtofile*namefile
            CSV.write(fullpath, df)
        end
    end
end