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

using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

include("/home/sergey/work/repo/dynamical-systems/FHN_Korotkov/PDF_clear_version/detect_spike.jl")

function two_coupled_rulkov_first_iteration(u, p)

    function right_part_x(x, y, z)
        if x <= 0.0
            return α/(1.0 - x) + y
        elseif (0.0 < x < (α + y)) && (z <=0.0)
            return α  +y
        elseif (x >= α + y) || (z > 0)
            return -1.0
        else
            return -1.0;
        end
    end

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, x2, y2, z2 = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, k, γ_1, γ_2, g1, g2 = p

    I21 = g2 * ( x_rp - x1 ) * xi(x2)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21, z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21)
    z1n = x1

    I12 = g1 * ( x_rp - x2 ) * xi(x1)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12, z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12)
    z2n = x2

    return SVector(x1n, y1n, z1n, I21, x2n, y2n, z2n, I12)
end

function two_coupled_rulkov(u, p, t)

    function right_part_x(x, y, z)
        if x <= 0.0
            return α/(1.0 - x) + y
        elseif (0.0 < x < (α + y)) && (z <=0.0)
            return α  +y
        elseif (x >= α + y) || (z > 0)
            return -1.0
        else
            return -1.0;
        end
    end

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, I21prev, x2, y2, z2, I12prev = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, k, γ_1, γ_2, g1, g2 = p

    I21 = γ_2 * I21prev + g2 * ( x_rp - x1 ) * xi(x2)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21, z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21)
    z1n = x1

    I12 = γ_1 * I12prev + g1 * ( x_rp - x2 ) * xi(x1)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12, z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12)
    z2n = x2

    return SVector{8}(x1n, y1n, z1n, I21, x2n, y2n, z2n, I12)
end

function get_params_two_coupled_rulkov()
    α = 4.6; σ =-0.1; μ = 0.001;
    β_syn = 0.0001; σ_syn = 1.0; k = 1.0;
    x_rp = -1.5; x_th = -0.8;
    γ_1 = 0.0; γ_2 = 0.0;
    g1 = 0.8; g2 = 0.2;
    return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, k, γ_1, γ_2, g1, g2]
end

function pdf(array_peaks, threshold_range, ϵ)
    total_number_peaks = length(array_peaks)
    array_PDF = Float64[]
    for index in range(1, length(threshold_range), step = 1)
        number_peaks_crossing_theshold = length(findall( (threshold_range[index] + ϵ) .>= array_peaks .>= threshold_range[index] ) )
        PDF_for_theshold = number_peaks_crossing_theshold / total_number_peaks
        push!(array_PDF, PDF_for_theshold)
    end
    return array_PDF
end

function get_PDF(X, shift)
    total_count_IEI = length(X)
    array_PDF = Float64[]
    for index in 1:total_count_IEI
        count_IEI_i = count(X[index]-shift .<= X .<= X[index]+shift)
        PDF_IEI_i = count_IEI_i / total_count_IEI
        push!(array_PDF, PDF_IEI_i)
    end
    return array_PDF
end

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

params = get_params_two_coupled_rulkov();
params[11] = 0.8;
params[12] = 0.2;
u0 = [1.0, 0.3, 0.01, 2.0, 0.7, 0.6];
tspan = (0, 1000000);
u0_first_iteration = two_coupled_rulkov_first_iteration(u0, params);

prob = DiscreteProblem(two_coupled_rulkov, SVector{8}(u0_first_iteration), tspan, params);
sol = solve(prob);
t_start = 200000;
x_sum = sol[1, t_start:end] + sol[5, t_start:end]
t_range = sol.t[t_start:end]

array_local_max, array_t_local_max = get_local_max([x_sum, t_range])

index_local_max_above_zero = findall(x-> x >= 0.0, array_local_max)
peaks_local_max_above_zero = array_local_max[index_local_max_above_zero]

Hs_x_sum = Hs(array_local_max, 6)

threshold_range = range(0, 3, length = 1000000);
ϵ = 0.1;
array_pdf = pdf(array_local_max, threshold_range, 1.0);

f = Figure()
ax = Axis(f[1, 1], xlabel = L"peaks", ylabel = L"PDF_{peaks} OLD", yscale = log10)
lines!(ax, threshold_range, array_pdf, linewidth = 1.0, color = :blue)
vlines!(ax, Hs_x_sum, linewidth = 3.0, linestyle = :dash, color = :red)
display(GLMakie.Screen(), f)


t_plot_end = 10000 #length(t_range) #1_000_000;

f = Figure(size = (1000, 400))
ax = Axis(f[1, 1])
lines!(ax, t_range[1:t_plot_end], x_sum[1:t_plot_end], linewidth = 1.0, color = :black)
hlines!(ax, Hs_x_sum, linewidth = 3.0, linestyle = :dash, color = :red)
scatter!(ax, array_t_local_max, array_local_max, markersize = 7.0, color = :orange)
xlims!(ax, t_range[1], t_range[t_plot_end])
display(GLMakie.Screen(), f)


#= array_local_max_sort = sort(peaks_local_max_above_zero)
array_pdf_v2 = get_PDF(array_local_max_sort, ϵ)

f = Figure()
ax = Axis(f[1, 1], xlabel = L"peaks", ylabel = L"PDF_{peaks}", yscale = log10)
lines!(ax, array_local_max_sort, array_pdf_v2, linewidth = 1.0, color = :blue)
vlines!(ax, Hs_x_sum, linewidth = 3.0, linestyle = :dash, color = :red)
display(GLMakie.Screen(), f) =#

