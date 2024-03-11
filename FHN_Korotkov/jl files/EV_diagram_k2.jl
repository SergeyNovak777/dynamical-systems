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

using DifferentialEquations, DynamicalSystems, StaticArrays, JLD2

function EV_diagram(path_to_u0s, name,
    sys, params, index_control_parameter, range_control_param,
    )
    array_u0 = load(path_to_u0s)[name]


    
end

#= 
path = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/u0s_k2_length_2000_for_direction_1.jld2"
name = "u0ss"

function peaks(x)
    peaks_ = Float64[]
    len_ = length(x)
    for i in range(2, len_ - 1, step = 1)
        if x[i-1] < x[i] > x[i+1]
            push!(peaks_, x[i])
        end
    end
    return peaks_
end

function calc_number_EEs(x)
    threshold = Hs(x)
    counts = length(x[x.>=threshold])
    return counts
end

Hs_ = Statistics.mean(xsum) + 6 * Statistics.std(xsum)
=#