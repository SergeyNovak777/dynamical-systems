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
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/pdf_function.jl")
using StaticArrays, DifferentialEquations, DynamicalSystems, Statistics, CairoMakie

function symmetrical(u, p, t)
    x, y = u
    f,ω_e, α, λ, ω_o    = p

    xdt = y
    ydt = (f * cos(ω_e * t) - α * y - λ * x * y^2 - ω_o * x) / (1 + λ * x^2)
    return SVector(xdt, ydt)
end

function get_params_mems()
    f = 3.1665; ω_e = 1.0; α = 0.2; λ = 0.5; ω_o = 0.25;
    return [f,ω_e, α, λ, ω_o]
end

u0 = [2.9347001488524063, -0.3563306330533323]
tspan = (0.0, 300000.0)
parameters = get_params_mems()
prob = ODEProblem(symmetrical, SVector{length(u0)}(u0), tspan, parameters)
sol = solve(prob, RK4(), adaptive = true, abstol = 1e-9, reltol = 1e-9, maxiters = 1e8)
len_sol = length(sol)

data = [sol[1, :], sol.t]
array_peaks, array_t_peaks = get_peaks(data; level_zero = "above")
Hs_x = Hs_above(array_peaks, 4)
σ = Statistics.mean(array_peaks)

width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
lines!(ax, sol.t, sol[1, :], color = :blue, linewidth = 0.5)
hlines!(Hs_x, color = "red", linewidth = 2.0, linestyle = :dash)
hlines!(σ , color = "lime", linewidth = 2.0, linestyle = :dash)
display(f)

count_thesholds = 1000000
maxvalue = maximum(sol[1, :])
threshold_range, array_PDF = pdf_v1(array_peaks, count_thesholds, maxvalue)

f = Figure(size = (400, 400))
ax = Axis(f[1, 1], yscale = log10)
lines!(threshold_range, array_PDF, linewidth = 3.0)
vlines!(ax, Hs_x, color = "red", linestyle = :dash, linewidth = 3.0)
display(f)

sel_1, pdf_1 = CALCPDF(array_peaks, count_thesholds)

f = Figure(size = (400, 400))
ax = Axis(f[1, 1], yscale = log10)
lines!(sel_1, pdf_1, linewidth = 3.0)
vlines!(ax, Hs_x, color = "red", linestyle = :dash, linewidth = 3.0)
display(f)

#=
Hs_above(x, k) = Statistics.mean(x) + k * Statistics.std(x)
Hs_below(x, k) = Statistics.mean(x) - k * Statistics.std(x)

function get_peaks(data; level_zero = "above")
    xrange, trange = data
    if level_zero == "above"
        xrange, trange = xrange[findall(xrange.>=0)], trange[findall(xrange.>=0)]
    else
        xrange, trange = xrange[findall(xrange.<=0)], trange[findall(xrange.<=0)]
    end
    array_peaks = Float64[]
    array_t_peaks = Float64[]
    for index in range(2, length(xrange)-1 )
        if xrange[index-1] < xrange[index] > xrange[index+1]
            push!(array_peaks, xrange[index])
            push!(array_t_peaks, trange[index])
        end
    end
    return array_peaks, array_t_peaks
end

threshold_range = range(0.0, maximum(sol[1, :]), length = count_thesholds)
total_number_peaks = length(array_peaks)
array_PDF = Float64[]
ϵ = 0.9
for index in range(1, count_thesholds, step = 1)
    number_peaks_crossing_theshold = length(findall( (threshold_range[index] + ϵ) .>= array_peaks .>= threshold_range[index] ) )
    PDF_for_theshold = number_peaks_crossing_theshold / total_number_peaks
    push!(array_PDF, PDF_for_theshold)
end

function CALCPDF(spikes, count_thesholds)
    thesholds = range(0,maximum(spikes),count_thesholds)
    ee_counter = [sum(i-> s<=i<s+0.9, spikes) for s in thesholds]
    pdf = ee_counter ./ length(spikes)
    return thesholds, pdf
end
=#