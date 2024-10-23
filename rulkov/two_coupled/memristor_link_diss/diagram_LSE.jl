username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

include("/home/sergey/work/repo/dynamical-systems/system.jl");
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl");
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl");

t_tr = (0, 3_000_000);
t_calc_LSE = 1_000_000;

u0 = SVector([-2.083638390440308, -3.9302148554862937, -2.0867818075436624,
            -1.97137793066347, -3.819163877171352, -1.9745518175123469,
            -0.11222999003131551])

prob = DiscreteProblem(rulkov_two_coupled_chem_mem, SVector{7}(u0), tspan, params);

# first iteration
sol = 

length_range_g2 = 500;
range_g2 = range(0.0, 10.0, length = length_range_g2);
array_LSEs = zeros(length_range_g2, length(u0));
array_u0s = zeros(length_range_g2, length(u0)); 
