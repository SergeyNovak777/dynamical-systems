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

using DifferentialEquations, DynamicalSystems, StaticArrays, GLMakie, JLD2

function FHN2_try3(u, p ,t)
    x1, y1, x2, y2, z= u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))
    ρz = k1 + k2 * z ^ 2

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + ρz * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + ρz * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    dzdt = x1 - x2
    return SVector(dx1dt, dy1dt, dx2dt, dy2dt, dzdt)
end

function FHN2_try3_params()
    ϵ = 0.01; a = -1.01;
    g = 0.1; k = 50.0; σ = 50.0 * pi / 180; α = 160.0 * pi / 180;
    k1 = 0.0; k2 = 0.0
    return [ ϵ, a, g, k, σ, α, k1, k2]
end

function plot_bifurcation_diagram(output, rangek1)
    markersize = 0.5;
    lbsize = 30;
    ticksize = 30;
    fig = Figure(size = (1200, 350))
    axis = Axis(fig[1,1],
    xlabel = L"k_1",  ylabel = L"x_1",
    xlabelsize = lbsize, ylabelsize = lbsize,
    xticklabelsize = ticksize,yticklabelsize = ticksize)
    
    for (j, p) in enumerate(rangek1)
        scatter!(axis, fill(p, length(output[j])), output[j];
         color = ("black", 0.5), markersize = markersize)
    end
    display(GLMakie.Screen(), fig)
end
function bifurcation_diagram()

    #u0 = [1.0, 0.0, 0.01, -1.0, 0.0]
    u0 = [-1.0, 0.0, 0.01, -1.0, 0.0]
    params = FHN2_try3_params()
    integ_set = (alg = RK4(), adaptive = false, dt = 0.001)

    ds = CoupledODEs(FHN2_try3, u0, params, diffeq = integ_set)
   
    t = 1000
    ttr = 500

    k1_start = 0.0
    k1_end = 0.094
    len = 3000
    rangek1 = range(k1_start, k1_end, length = len)
    index_control_param = 7

    index_saving_var = 1
    index_fixed_var = 4
    value_fixed_var = 0.0
    surface = (index_fixed_var, value_fixed_var)
    setting_root = (xrtol = 1e-11, atol = 1e-11)
    pmap = PoincareMap(ds, surface, rootkw = setting_root)
    
    output = orbitdiagram(pmap, index_saving_var, index_control_param, rangek1;
     n = t, Ttr = ttr, show_progress = true)

     return output
end

k1_start = 0.0
k1_end = 0.094
len = 3000
rangek1 = range(k1_start, k1_end, length = len)

output = bifurcation_diagram()

plot_bifurcation_diagram(output, rangek1)

jld2save("bif_dia_k1_length_3000.jld2";output)