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

using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

function FHN2_4d(u, p ,t)
    x1, y1, x2, y2 = u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + (k1 + k2 * (y1 - y2)^2) * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + (k1 + k2 * (y1 - y2)^2) * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    return SVector(dx1dt, dy1dt, dx2dt, dy2dt)
end

u0 = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];
params = FHN2_try3_params()
params[8] = 64.76190476190476; 
integ_set = (alg = Vern9(), adaptive = true, abstol=1e-13, reltol=1e-13, maxiters = 1e8)

ds = CoupledODEs(FHN2_4d, u0, params, diffeq = integ_set)

t = 2000
ttr = 1000

k1_start = 0.09686
k1_end = 0.09
len = 200
rangek1 = range(k1_start, k1_end, length = len)
index_control_param = 7

index_saving_var = 1
index_fixed_var = 3
value_fixed_var = -1.01
surface = (index_fixed_var, value_fixed_var)
setting_root = (xrtol = 1e-11, atol = 1e-11)
pmap = PoincareMap(ds, surface, rootkw = setting_root)

output = orbitdiagram(pmap, index_saving_var, index_control_param, rangek1;
 n = t, Ttr = ttr, show_progress = true)


 markersize = 0.5;
lbsize = 30;
ticksize = 30;
fig = Figure(size = (1200, 350))
axis = Axis(fig[1,1],
xlabel = L"k_1",  ylabel = L"x_1",
xlabelsize = lbsize, ylabelsize = lbsize,
xticklabelsize = ticksize,yticklabelsize = ticksize)

for (j, p) in enumerate(rangek1)
scatter!(axis, fill(p, length(output[j])), output[j]; color = ("black", 0.5), markersize = markersize)
end
display(GLMakie.Screen(), fig)