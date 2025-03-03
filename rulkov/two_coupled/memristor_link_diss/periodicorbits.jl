username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

include("/home/sergey/work/repo/dynamical-systems/system.jl");

params = get_params_rulkov_two_coupled_chem_mem()

params[1] = 3.9; # α
params[2] = 1.0; # σ
params[10] = 5.0; # g1
params[11] = 1.0; # g2

params[12] = 0.1; # k1
params[13] = 0.00; # k2

u0 = SVector{7}(-1.0, -2.5729418525243473, 1.1604540938548373, -0.8654368209270067, -3.050292978611296, -0.9789719381322224, 2.1394260319870595);

ds = DeterministicIteratedMap(rulkov_two_coupled_chem_mem, SVector{7}(u0), params)

xs = range(-5, 5, length = 4);
ys = range(-3, 3, length = 4);
zs = range(-5, 5, length = 4);
ls = range(-3, 3, length = 4);

ics = [SVector{7}(x1, y1, z1, x2, y2, z2 , l) for x1 in xs for y1 in ys for z1 in zs for x2 in xs for y2 in ys for z2 in zs for l in ls]

periodicorbits(ds, 6, ics)