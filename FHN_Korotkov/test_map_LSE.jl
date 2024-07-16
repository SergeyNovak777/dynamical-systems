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

include("/home/sergey/work/repo/dynamical-systems/Map_and_diagram_LSE/map_LSE_linear.jl");

using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

sys = FHN2_try3;
params = parameters = FHN2_try3_params();
u0 = [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35];
u0 = SVector{5}(u0);

alg = Vern9();
adaptive = true;
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

Ttr = 3_000;
t_LSE = 10_000;

len = 10;
name_p1 = "k1";
name_p2 = "k2";
index_p1 = 7;
index_p2 = 8;
range_p1 = range(0.0, 0.08, length = len);
range_p2 = range(0.0, 10.0, length = len);

type_inheritance = "from_left_to_right";

map_LSE_linear(sys, parameters, u0, integrator_setting,
    Ttr, t_LSE,
    name_p1, name_p2, index_p1, index_p2, range_p1, range_p2,
    type_inheritance)