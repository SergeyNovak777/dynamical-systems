user = "Alex"
pathtorepo = "C:\\Users\\" * user * "\\Desktop\\"
using Pkg
Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
include("C:\\Users\\" * user * "\\Desktop\\dynamical-systems\\system.jl")

using StaticArrays, DifferentialEquations, DynamicalSystems, Statistics, CairoMakie

# setting method of integrate
tstep = 0.01
integ_set = (alg = RK4(), adaptive = false, dt = tstep)

# parameters
a = 1.0; b = 3.0; c = 1.0; d = 5.0;
xr = -1.6; r = 0.01; s = 5.0; I = 4.0; xv = 2.0;
k1= -0.17; k2 = -0.17;
k1_me = 0.0
k2_me = 0.0

u0 = [-1.5, 0.0, 0.0, -2.5, 0.0, 0.0, 0.0]
p = [a, b, c, d, s, xr, r, I, xv, k1, k2, k1_me, k2_me];
ds = CoupledODEs(HR_mem, u0, p, diffeq = integ_set)

k2_merange = range( 0.0, 1.0, length = 100 )

idx_control_parameter = 12
idx_save = 1
idx_fix = 4; fixed_value = 0.0
surface = (idx_fix, fixed_value)
setting_root = (xrtol = 1e-15, atol = 1e-20)


pmap = PoincareMap(ds, surface, rootkw = setting_root)
output = orbitdiagram(pmap, idx_save, idx_control_parameter, k2_merange; n = 500, Ttr = 500, show_progress = true)