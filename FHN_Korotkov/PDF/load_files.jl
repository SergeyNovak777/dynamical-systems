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
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/pdf_function.jl")
    include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/preprocessing_spike.jl")
end

using StaticArrays, Statistics, CairoMakie, GLMakie, JLD2

path_to_save = "/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF/sol_k2=75.74/"

EE_peaks = Float64[]
EE_times = Float64[]
sol = Float64[]
sol_times = Float64[]

function vcat_arrays(sol, sol_times, EE_peaks, EE_times)
    for iteration in range(1, 10)

        namefile_sol_x1 = "$(iteration)_sol_x1.jld"
        namefile_EE = "$(iteration)_EE.jld"
        full_path_to_save_sol_x1 = path_to_save*namefile_sol_x1
        full_path_to_save_EE = path_to_save*namefile_EE

        sol_x1, sol_t = load(full_path_to_save_sol_x1)["data_x1"]
        EE_peaks_and_times = load(full_path_to_save_EE)
        peaks_EE = EE_peaks_and_times["peaks_EE"]
        t_EE = EE_peaks_and_times["t_EE"]

        EE_peaks = vcat(EE_peaks, peaks_EE)
        EE_times = vcat(EE_times, t_EE)
        sol = vcat(sol, sol_x1)
        sol_times = vcat(sol_times, sol_t)
    end

    return sol, sol_times, EE_peaks, EE_times
end

sol, sol_times, EE_peaks, EE_times = vcat_arrays(sol, sol_times, EE_peaks, EE_times)
length(sol)

tstart = 1; tend = Int64(6e6)
width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
lines!(ax, sol_times[tstart:tend], sol[tstart:tend])
scatter!(ax, EE_times, EE_peaks, markersize = 10, color = :red)
xlims!(sol_times[tstart], sol_times[tend])
display(GLMakie.Screen(), f)