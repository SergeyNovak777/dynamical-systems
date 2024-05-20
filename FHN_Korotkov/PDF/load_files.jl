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

path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/timeseries_k2_75/"

len_part_timeseris = Int64[]
sol = Float64[]
sol_times = Float64[]

function vcat_arrays(sol, sol_times)
    for iteration in range(1, 27)

        namefile_sol_x1 = "$(iteration)_sol_x1.jld"
        full_path_to_save_sol_x1 = path_to_save*namefile_sol_x1

        sol_x1, sol_t = load(full_path_to_save_sol_x1)["data_x1"]
        sol = vcat(sol, sol_x1)
        sol_times = vcat(sol_times, sol_t)
        push!(len_part_timeseris, length(sol_x1))
    end

    return sol, sol_times
end

sol, sol_times = vcat_arrays(sol, sol_times)
println(length(sol))
#-------------------------------------------------------------

data_x1 = [sol, sol_times]
sol = nothing; sol_times = nothing

array_spikes_max, array_t_spikes_max = get_peaks(data_x1; level_zero = "nothing")
array_spikes_thresholds, array_t_spikes_thresholds = get_peaks_neg(data_x1)
baseline = Statistics.mean(array_spikes_thresholds)
println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")

if check_timeseries(array_spikes_max, array_spikes_thresholds) == true
    array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
else
    drop_false_start_end(array_t_spikes_max, array_spikes_max, array_t_spikes_thresholds)
    array_t_spikes_max_correct, array_spikes_max_correct, amplitudes = select_spike_by_amplitude(array_spikes_thresholds, array_t_spikes_max, array_spikes_max)
end

println("length local max: $(length(array_spikes_max)); length local min $(length(array_spikes_thresholds))")
println("count correct spike: $(length(array_spikes_max_correct))")

Hs_x  = Hs_above(amplitudes, 8)
index_EE = findall(x-> x >= Hs_x, array_spikes_max_correct)
t_EE = array_t_spikes_max_correct[index_EE]
EE = array_spikes_max_correct[index_EE]
len_EE = length(EE)
println("count EE: $(len_EE)")

#----------------------------------------------
IEI = Float64[]
IEI_v2 = zeros(len_EE-1, 3)
for index in range(1, len_EE-1)
    t_EE_i = t_EE[index]
    t_EE_i_plus_one = t_EE[index+1]
    IEI_i = t_EE_i_plus_one - t_EE_i
    push!(IEI, IEI_i)
    IEI_v2[index, :] = [index, index+1, IEI_i]
end
IEI_v2 = IEI_v2[sortperm(IEI_v2[:, 3]), :]
IEI = sort(IEI)
IEI = IEI[IEI.>=10]
len_IEI = length(IEI)
ϵ = 10.0;

count_identical_IEI = zeros(len_IEI)
PDF_IEI = zeros(len_IEI)

for index in range(1, len_IEI)
    count_IEI_i = count(IEI[index] .<= IEI .<= IEI[index]+ϵ)
    count_identical_IEI[index] = count_IEI_i
    PDF_IEI_i = count_IEI_i / len_IEI
    PDF_IEI[index] = PDF_IEI_i
end
#--------------------------------------------
height_window = 400; width_window = 1100;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"IEI", ylabel =  L"PDF")
lines!(ax, IEI, PDF_IEI, color = :black, linewidth = 1.0)
display(GLMakie.Screen(), f)
#--------------------------------------------
#=
part_time_series = 50
tstart = 1; tend = Int64(5e6) #length(data_x1[2])
width_window = 500; height_window = 1000;
f = Figure(size = (height_window, width_window))
ax = Axis(f[1, 1], xlabel = L"time", ylabel =  L"x")
lines!(ax, data_x1[2][tstart:tend], data_x1[1][tstart:tend])
scatter!(ax, array_t_spikes_max_correct, array_spikes_max_correct, markersize = 10, color = :green)
#scatter!(ax, array_t_spikes_max, array_spikes_max, markersize = 10, color = :green)
scatter!(ax, t_EE, EE, markersize = 10, color = :deeppink)
scatter!(ax, array_t_spikes_thresholds, array_spikes_thresholds, markersize = 10, color = :blue)
xlims!(data_x1[2][tstart], data_x1[2][tend])
display(GLMakie.Screen(), f)=#