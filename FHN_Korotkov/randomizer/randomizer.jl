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
#INCLUDES
include("/home/sergey/work/repo/dynamical-systems/system.jl")
#-----------------------------------------------------------------
#PACKAGES
using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2
using DataFrames, CSV
#-----------------------------------------------------------------
function generate_u0(limits_X, limits_Y, len_array_u0)

    x1_2_range =range(limits_X[1], limits_X[2], length = len_array_u0)
    y1_2_range = range(limits_Y[1], limits_Y[2], length = len_array_u0)

    array_u0 = zeros(len_array_u0, 5)

    for index in range(1, len_array_u0)
        u0_x1_x2 = rand(x1_2_range, (2,1))
        u0_y1_y2 = rand(y1_2_range, (2,1))
        array_u0[index, 1], array_u0[index, 3] = u0_x1_x2
        array_u0[index, 2], array_u0[index, 4] = u0_y1_y2
        array_u0[index, 5] = u0_y1_y2[1] - u0_y1_y2[2]
    end

    return array_u0
end

function solver(sys, u0, params, integrator_setting, tspan)

    prob = ODEProblem(sys, SVector{length(u0)}(u0), tspan, params)

    if integrator_setting.adaptive == true
        sol = solve(prob, integrator_setting.alg, adaptive = true, abstol = integrator_setting.abstol, reltol  = integrator_setting.reltol)
    else
        sol = solve(prob, integrator_setting.alg, adaptive = false, dt = integrator_setting.dt, maxiters = integrator_setting.maxiters)
    end    
    return sol
end

function calculate_LSE(sys, u0, params, integrator_setting, time_calculate)
    ds = CoupledODEs(sys, u0, params,
    diffeq = integrator_setting);
    LSE = lyapunovspectrum(ds, time_calculate)
    return LSE
end

function create_df(u0, last_point, LSE)
    df = DataFrame("start x_1" => u0[1], "start y_1" => u0[2], " startx_2" => u0[3],
    "start y_2" => u0[4],"start z" => u0[5],
    "end x_1" => last_point[1], "end y_1" => last_point[2], " end x_2" => last_point[3],
    "end y_2" => last_point[4],"end z" => last_point[5],
    "LSE1" => LSE[1], "LSE2" => LSE[2], "LSE3" => LSE[3], "LSE4" => LSE[4], "LSE5" => LSE[5])
    return df
end

function push_df(df, u0, last_point, LSE)
    push!(df, (u0[1], u0[2], u0[3], u0[4], u0[5], last_point[1], last_point[2], last_point[3], last_point[4], last_point[5], LSE[1], LSE[2], LSE[3], LSE[4], LSE[5]))
end

function randomizer_u0(sys, params, tspan_solution, time_calculate_LSE,
    integrator_setting, len_matrix_random_u0, limits_X, limits_Y, path_to_save)

    matrix_random_u0 = generate_u0(limits_X, limits_Y, len_matrix_random_u0)

    for index in range(1, len_matrix_random_u0, step = 1)
        u0 = matrix_random_u0[index, :]
        solution = solver(sys, u0, params, integrator_setting, tspan_solution)
        last_point = solution[end]
        LSE = calculate_LSE(sys, last_point, params, integrator_setting, time_calculate_LSE)

        if index === 1
            global df = create_df(u0, last_point, LSE)
            CSV.write(path_to_save, df)
        else
            push_df(df, u0, last_point, LSE)
        end

        if mod(index, 10) == 0
            CSV.write(path_to_save, df)
        end
    end
end

function main()

    sys = FHN2_try3
    params = FHN2_try3_params()
    tspan_solution = 10000
    time_calculate_LSE = 20000
    params[3] = 0.02
    params[7] = 0.150
    params[8] = 0.0
    integrator_setting = (alg = Vern9(), adaptive = true,
    abstol = 1e-11, reltol = 1e-11)
    len_matrix_random_u0 = 1000
    limits_X = [-2.0, 2.0]
    limits_Y = [-2.0, 2.0]

    path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/randomizer/"
    file_name = "$len_matrix_random_u0" * "_random_u0_g_0_02_k1_$(params[7])_k2_$(params[8]).csv"
    fullpath = path_to_save * file_name

    randomizer_u0(sys, params, tspan_solution, time_calculate_LSE,
    integrator_setting, len_matrix_random_u0, limits_X, limits_Y, fullpath)

end