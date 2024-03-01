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

include("/home/sergey/work/repo/dynamical-systems/system.jl")

using DifferentialEquations, DynamicalSystems, StaticArrays, JLD2, GLMakie

function goto_attractor(u0_loc, parameters, tspan)
    prob = ODEProblem(FHN2_try3, u0_loc, tspan, parameters)
    sol = solve(prob, RK4(), adaptive = false, dt = 0.001, maxiters = 5e6,
    save_everystep = false, save_start = false);
    last_point = sol[end]
    return last_point
end

function save_files(filenameu0s, filenameLSEs)
    jldsave(filenameu0s; u0ss)
    jldsave(filenameLSEs; LSEs)
end

function calculate_LSE(u0_lc, params, time_calculate_LLE)
    ds = CoupledODEs(FHN2_try3, u0_lc, params,
     diffeq = ( alg = RK4(), adaptive = false, dt = 0.001, maxiters = 5e6 ));
    LLE = lyapunovspectrum(ds, time_calculate_LLE)
    return LLE
end

function print_output(index, parameter, u0, point, LSE)
    println("index: $(index); k_1: $(parameter)"); flush(stdout)
    println("u0: $(u0)"); flush(stdout)
    println("last point: $(point)");flush(stdout)
    println("LSE: $(LSE)"); flush(stdout)
    println("-----------------------------"); flush(stdout)
    println("");flush(stdout)
end

function inheritance_limit_cycle(index_p, range_param, tspan, time_calculate_LLE)

    u0 = SVector{5}([-1.0, 0.0, 0.01, -1.0, 0.0])
    parameters = FHN2_try3_params()
    
    for index in range(1,length(range_param))

        if index == 1
            global u0_local = u0
        end

        parameters[index_p] = range_param[index]
        point = goto_attractor(u0_local, parameters, tspan)
        LSE = calculate_LSE(point, parameters, time_calculate_LLE)

        print_output(index, parameters[index_p], u0_local, point, LSE)

        u0_local = point

        u0ss[index, :] = point
        LSEs[index, :] = LSE


        if mod(index, 10) == 0
            save_files(filenameu0s, filenameLSEs)
        end
    end
end

cd("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/data")

tspan = (0.0, 500)
time_calculate_LLE = 5000

k1_start = 0.0
k1_end = 0.094
len = 3000
index_param = 7
range_param = range(k1_start, k1_end, length = len)

u0ss = zeros(len, 5)
LSEs = zeros(len, 5)

filenameu0s = "u0s_k1_length_3000.jld2"
filenameLSEs = "LSEs_k1_length_3000.jld2"

inheritance_limit_cycle(index_param, range_param, tspan, time_calculate_LLE)

