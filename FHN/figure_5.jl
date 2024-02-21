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

function plot_bifurcation_diagram(output, rangek1)
    markersize = 5.0;
    lbsize = 40;
    ticksize = 35;
    fig = Figure(resolution = (1200, 700))
    axis = Axis(fig[1,1],
    xlabel = L"k_1",  ylabel = L"x_1",
    xlabelsize = lbsize, ylabelsize = lbsize,
    xticklabelsize = ticksize,yticklabelsize = ticksize)
    
    for (j, p) in enumerate(rangek1)
        scatter!(axis, fill(p, length(output[j])), output[j];
         color = ("blue", 1.0), markersize = markersize)
    end
    display(fig)
end
function bifurcation_diagram()

    u0 = [-1.0690237112785876, -0.5773325275931365, -1.1977568518748323, -0.5654141361024896, 1.60676442231909478]
    params = two_coupled_fhn_get_params()
    integ_set = (alg = RK4(), adaptive = false, dt = 0.001)

    ds = CoupledODEs(two_coupled_fhn, u0, params, diffeq = integ_set)
   
    t = 1000
    ttr = 500

    k1_start = 0.0
    k2_end = 0.08
    len = 100
    rangek1 = range(k1_start, k2_end, length = len)
    index_control_param = 7

    index_saving_var = 1
    index_fixed_var = 4
    value_fixed_var = 0.0
    surface = (index_fixed_var, value_fixed_var)
    pmap = PoincareMap(ds, surface)
    
    setting_root = (xrtol = 1e-11, atol = 1e-11)

    output = orbitdiagram(pmap, index_saving_var, index_control_param, rangek1;
     n = t, Ttr = ttr, show_progress = true)

     return output
end

k1_start = 0.0
k2_end = 0.08
len = 100
rangek1 = range(k1_start, k2_end, length = len)

output = bifurcation_diagram()

plot_bifurcation_diagram(output, rangek1)