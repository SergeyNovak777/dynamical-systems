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

function generate_u0(parameters, len_array_u0)

    range_s_E = range(parameters[3], 1.0, length = 100)
    range_s_I = range(parameters[6], 1.0, length = 100)
    range_r = range(0.0, 1.0, length = 100)
    range_astro = range(0.0, parameters[20] * parameters[19], length = 100)

    
    array_u0 = zeros(len_array_u0, 5)

    for index in range(1, len_array_u0)

        u0_1 = (rand(range_s_E, 1))[1]
        u0_2 = (rand(range_s_I, 1))[1]
        u0_3_4 = rand(range_r, (2,1))
        u0_5 = (rand(range_astro, 1))[1]

        array_u0[index, 1] = u0_1
        array_u0[index, 2] = u0_2
        array_u0[index, 3:4] = u0_3_4
        array_u0[index, 5] = u0_5

    end
    return array_u0
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
    len_array_u0 = 250
    array_u0 = generate_u0(parameters, len_array_u0)

    # index gammaY 21
    len_range_param = 1000
    range_gammaY = range(0.0, 100.0, length = len_range_param)

    for index_param_in_array in range(1, len_range_param, step = 1)

        parameters[21] = range_gammaY[index_param_in_array]

        for index in range(1, len_array_u0, step = 1)
            
            prob = init_ODE_prob(system, parameters, array_u0[index, :], t_attract)
            last_point = goto_attractor(prob, integrator_setting)
            ds = init_Coupled_ODE(system, parameters, last_point, integrator_setting)
            LSE = calculate_LSE(ds, t_LSE)

            sE_sv, sI_sv, rE_sv, rI_sv, Y_sv = array_u0[index, :]
            sE_ev, sI_ev, rE_ev, rI_ev, Y_ev = last_point
            l1, l2, l3, l4, l5 = LSE
            
            if index == 1 && index_param_in_array == 1
                global df = DataFrame("gammaY" => parameters[21], "sE_start" => sE_sv, "sI_start" => sI_sv, "rE_start" => rE_sv, "rI_start" => rI_sv, "Y_start" => Y_sv,
                    "sE_end" => sE_ev, "sI_end" => sI_ev, "rE_end" => rE_ev, "rI_end" => rI_ev, "Y_end" => Y_ev,
                    "LE1" => l1, "LE2" => l2, "LE3" => l3, "LE4" => l4, "LE5" => l5)
            else
                push!(df, (parameters[21], sE_sv, sI_sv, rE_sv, rI_sv, Y_sv, sE_ev, sI_ev, rE_ev, rI_ev, Y_ev, l1, l2, l3, l4, l5))
            end

            if mod(index, 10) == 0
                pathtofile = "/home/sergey/work/repo/dynamical-systems/rate model/randomizer/"
                namefile = "gammaY_1000_u0s_250.csv"
                fullpath = pathtofile*namefile
                CSV.write(fullpath, df)
            end
        end
    end
         
end