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
end

using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

function FHN2_4d(u, p ,t)
    x1, y1, x2, y2 = u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + (k1 + k2 * (y1 - y2)^2) * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + (k1 + k2 * (y1 - y2)^2) * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    return SVector(dx1dt, dy1dt, dx2dt, dy2dt)
end

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9();
adaptive = true;
abs_tol = 1e-13;
rel_tol = 1e-13;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[3] = 0.1;
parameters[7] = 0.09353383 #0.09686; #0.094589;
parameters[8] = 64.7619047619048;

u0_start = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];
u0_start = SVector{4}(u0_start);

t_end = 10_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, sol[end], tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_4d, sol[end], parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 50000);
println(LSE);

labelsize = 50;
ticksize = 25;

length_sol = length(sol);
ttr = t_truncate(length_sol)
t_plot_start =  ttr
t_plot_end = t_plot_start + 10_000;
                
indexx = 1; indexy  = 3; indexz = 2;
f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_1",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false)#,
lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);


#= x_mean = (sol[1, ttr:end] + sol[3, ttr:end]) / 2;
t_start_plot_timeseries = 1; t_end_plot_timeseries = 500_000;
labelsize = 85;
ticksize = 50;

t_start_plot_timeseries = 1; t_end_plot_timeseries = 400_000;
f = Figure(size = (1200, 700));
ax_x_mean = Axis(f[1, 1], ylabel = L"\overline{x}",
    xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, xticklabelsvisible = false,
    xticks = [0, 2000, 4000])

ax_x1_2 = Axis(f[2, 1], xlabel = L"t", ylabel = L"x_i",
xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize,
xgridvisible = false, ygridvisible = false)

lines!(ax_x_mean, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], x_mean[t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :blue)

lines!(ax_x1_2, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], sol[1, t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :red)
lines!(ax_x1_2, sol.t[t_start_plot_timeseries:t_end_plot_timeseries], sol[3, t_start_plot_timeseries:t_end_plot_timeseries],
linewidth = 1.0, color = :green)
display(f); =#