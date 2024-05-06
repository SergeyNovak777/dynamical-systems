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
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/pdf_function.jl")
using StaticArrays, DifferentialEquations, DynamicalSystems, Statistics, CairoMakie, GLMakie

t_truncate(t) = floor(Int64, t / 2)

parameters = FHN2_try3_params()
tspan = (0.0, 50000.0)
parameters[7] = 0.09
parameters[8] = 75.74
u0 = [1.7, 0.7, -1.4, 0.35, 0.7 - 0.35]
u0 = SVector{5}(u0)
prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
sol = solve(prob, DP8(), adaptive = true, abstol = 1e-11, reltol = 1e-11, maxiters = 1e8);
len_sol = length(sol.t)
tstart = t_truncate(len_sol); tend = len_sol

data = [sol[1, tstart:tend], sol.t[tstart:tend]]
array_peaks, array_t_peaks = get_peaks(data; level_zero = "nothing")
# get_peaks_neg(data)
array_peaks_abs = abs.(array_peaks)
Hs_x = Hs_above(array_peaks_abs, 8)
σ = Statistics.mean(array_peaks_abs)

width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
lines!(ax, sol.t[2600000:tend], sol[1, 2600000:tend], color = :blue, linewidth = 0.5)
scatter!(ax, array_t_peaks, array_peaks, markersize = 10, color = :red)
xlims!(ax, 49700, 49930)
hlines!(Hs_x, color = "red", linewidth = 2.0, linestyle = :dash)
hlines!(σ , color = "lime", linewidth = 2.0, linestyle = :dash)
display(GLMakie.Screen(), f)    

count_thesholds = 1000000
maxvalue = maximum(array_peaks_abs) # minimum(sol[1, tstart:tend])
minvalue = minimum(array_peaks_abs)
ϵ = 1.5
threshold_range, array_PDF = pdf_v2(array_peaks_abs, count_thesholds, minvalue, maxvalue, ϵ)

f = Figure(size = (400, 400))
ax = Axis(f[1, 1])
lines!(threshold_range, array_PDF, linewidth = 1.5)
vlines!(ax, Hs_x, color = "red", linestyle = :dash, linewidth = 1.5)
display(GLMakie.Screen(), f)


f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
scatter!(ax, array_t_peaks, abs.(array_peaks), markersize = 10, color = :red)
xlims!(ax, 49000, 49930)
hlines!(abs.(Hs_x), color = "red", linewidth = 2.0, linestyle = :dash)
display(GLMakie.Screen(), f)   
