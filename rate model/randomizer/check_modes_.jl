if Sys.iswindows()
    cd("D:\\work\\dynamical-systems\\env\\")
    using Pkg
    Pkg.activate("integrate")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
end
#-----------------------------------------------------------------
#INCLUDES
include(pathtorepo*"/system.jl")
#-----------------------------------------------------------------
#PACKAGES
using DifferentialEquations, DynamicalSystems, CSV, DataFrames, DelimitedFiles
using GLMakie
#-----------------------------------------------------------------

function preproc_df()
    
    df = DataFrame(CSV.File("/home/sergey/work/repo/dynamical-systems/rate model/randomizer/50000_u0s.csv"))
    df_filter = filter(r -> r.LE1 >= 0, df)

    u0s = df_filter[:,[:"sE_start",:"sI_start", :"rE_start", :"rI_start", :"Y_start"]]

    u0s = Matrix{Float64}(u0s) 
    return u0s
end

function plot_phase_space(X, labels, window_size, lw)
    x1, x2, x3 = X
    x1label, x2label, x3label = labels

    f = Figure(size = (window_size[1], window_size[2]))
    ax = Axis3(f[1, 1], xlabel = x1label, ylabel = x2label, zlabel = x3label)
    lines!(ax, x1, x2, x3, linewidth = lw)
    display(GLMakie.Screen(), f)
end

function main(index)

    data_from_csv = preproc_df()
    u0 = data_from_csv[index, :]

    sys = rate_model
    params = rate_model_get_params()
    t_attract = 1000.0
    integrator_setting = (alg = Vern9(), adaptive = true, abstol = 1e-11, reltol = 1e-11)
    tspan = (0.0, t_attract)
   
    prob = ODEProblem(sys, SVector{length(u0)}(u0), tspan, params)
    sol = solve(prob, alg = integrator_setting.alg,
        adaptive = integrator_setting.adaptive,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol)

    t_start_plot = floor(Int64, length(sol[1, :] )/ 2); t_end_plot = length(sol[1, :])
    X = [sol[1, t_start_plot:t_end_plot], sol[2, t_start_plot:t_end_plot], sol[3, t_start_plot:t_end_plot]]
    labels = ["s_E", "s_I", "r_E"]
    window_size = [600, 400]
    lw = 1.0

    plot_phase_space(X, labels, window_size, lw)

end