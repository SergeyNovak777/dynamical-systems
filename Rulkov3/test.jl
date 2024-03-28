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
#---------------------------------------------------------------------
# PACKAGES
using DynamicalSystems, StaticArrays, CairoMakie
#---------------------------------------------------------------------
#INCLUDE
include("/home/sergey/work/repo/dynamical-systems/Rulkov3/system.jl")

sys = rulkov
u0 = [0.0, 0.0]
t = 10000
params = rulkov_get_params()
params[1] = 3.9
params[2] = 0.04

ds = DeterministicIteratedMap(sys, u0, params)

tr, trange = trajectory(ds, t)

tstart = 6000; tend = 10000

f = Figure(size = (1000, 250))
ax = Axis(f[1,1], xlabel = "time", ylabel = "x")
lines!(ax, trange[tstart:tend], tr[tstart:tend,1])
display(f)