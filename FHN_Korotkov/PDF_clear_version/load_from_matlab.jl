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

using Statistics, MAT, JLD2, Plots

include("/home/sergey/work/repo/dynamical-systems/system.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")
include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/IEI.jl")
#include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/plotter.jl")

t_truncate(t) = floor(Int64, t / 2)
Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function load_vcat(sol, sol_times, len_charts)
    for index in 1:len_charts
        
        println("index: $index"); flush(stdout)
        filename_sol = "$(index)_solution.mat"
        filename_sol_t = "$(index)_t_interval.mat"

        println("open files"); flush(stdout)
        solve = matopen(folder*filename_sol)
        solve_t = matopen(folder*filename_sol_t)
        solve = read(solve, "solve")
        solve_t = read(solve_t, "t_solve")

        
        if index == 1
            len_sol = length(solve_t)
            ttr = t_truncate(len_sol)
            solve = solve[ttr:len_sol]
            solve_t = solve_t[ttr:len_sol]
        end

        println("vcat"); flush(stdout)
        sol = vcat(sol, solve);
        sol_times = vcat(sol_times, solve_t);

        println("-----------------------"); flush(stdout)
    end

    return sol, sol_times
end

folder = "/home/sergey/MEGA/matlab/files_long_time_series/"

println("index: 30"); flush(stdout)
filename_sol = "30_solution.mat"
filename_sol_t = "30_t_interval.mat"

println("open files"); flush(stdout)
solve = matopen(folder*filename_sol)
solve_t = matopen(folder*filename_sol_t)
solve = read(solve, "solve")
solve_t = read(solve_t, "t_solve")

len_solve = length(solve_t)
tstart = len_solve - 10000; tend = len_solve
plot(solve_t[tstart:tend], solve[tstart:tend])


#= f = Figure()
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_1")
lines!(ax, solve_t[tstart:tend], solve[tstart:tend], linewidth = 1.0, color = :blue)
display(f) =#

#= sol = Float64[]
sol_times = Float64[]

len_charts = 30

sol, sol_times = load_vcat(sol, sol_times, len_charts)

data = [sol, sol_times];
sol = nothing;
sol_times = nothing;
GC.gc();

data_local_max = get_local_max(data)
data_local_min = get_local_min(data)
drop_artifacts(data_local_max, data_local_min)

all_amplitudes = get_amplitudes_all_events(data_local_max[1], data_local_min[1])
mean_amplitudes = Statistics.mean(all_amplitudes)

peaks_spikes, t_peaks_spikes, amplitudes_above_mean = select_spikes(data_local_min[1], data_local_max, mean_amplitudes) # Statistics.std(all_amplitudes))
println("count spikes: $(length(peaks_spikes))")
Hs_x = Hs(amplitudes_above_mean, 8)

index_EEs = findall(x-> x >= Hs_x, peaks_spikes)
peaks_EEs = peaks_spikes[index_EEs]
t_EEs = t_peaks_spikes[index_EEs]
println("count EEs: $(length(peaks_EEs))")

array_IEI = get_IEI(t_EEs)
array_IEI = sort(array_IEI)
array_PDF_IEI = get_PDF_IEI(array_IEI; shift = 20)

plot_PDF_hist(array_IEI, array_PDF_IEI, 25)

len_solve = length(data[2])
tstart = len_solve - 10000; tend = len_solve
f = Figure()
ax = Axis(f[1, 1], xlabel = L"time", ylabel = L"x_1")
lines!(ax, data[2][tstart:tend], data[1][tstart:tend], linewidth = 1.0, color = :blue)
display(f) =#