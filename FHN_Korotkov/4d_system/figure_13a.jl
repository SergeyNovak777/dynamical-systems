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

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9();
adaptive = true;
abs_tol = 1e-13;
rel_tol = 1e-13;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[3] = 0.1;
parameters[7] = 0.09; #0.09;
parameters[8] = 76;   ; # 75.7;

u0_start = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];
u0_start = SVector{4}(u0_start);

t_end = 50_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters)
FHN2_4d([-1.01, -0.6367552038497962, -1.01, -0.6367552038515], parameters, 0)


sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);



path_to_save = "/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/images/rewrite_images/"
filename_hist = "fig_13_a_phase_space.eps"

labelsize = 80;
ticksize = 50;
CairoMakie.activate!();
length_sol = length(sol);
ttr = t_truncate(length_sol)
t_plot_start =  ttr
t_plot_end = t_plot_start + 100_000;
                
indexx = 1; indexy  = 3; indexz = 4;
f = Figure(size = (1000 ,600));
ax = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_2",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false,
    xlabeloffset = 60, ylabeloffset = 60, zlabeloffset = 85, protrusions = (20, 20,100, 20))

lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);
#save(path_to_save*filename_hist, f)


ds = CoupledODEs(FHN2_4d, sol[end], parameters,
diffeq = integrator_setting);
sol = nothing; GC.gc();

LSE = lyapunovspectrum(ds, 100_000);
println("LSE: $(LSE)");


#= x_intervals = interval(-1.5, 0);
y_intervals = interval(-1.5, 1.5);
box = x_intervals × y_intervals × x_intervals × y_intervals;
print("start calc fp");
fp, eig, _ = fixedpoints(ds, box, tol = 1e-12); =#

pmap = PoincareMap(ds, (4, -0.625));#0.0))

tr, trange = trajectory(pmap, 400_000)

tstartpo = 100_000; tendpo= 400_000;

f = Figure(size = (1000, 600))
ax = Axis(f[1, 1], xgridvisible = false, ygridvisible = false,
xlabel = L"x_1", ylabel = L"x_2", xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize)
scatter!(tr[tstartpo:tendpo, 1], tr[tstartpo:tendpo, 3], color = :red, markersize = 1.0)
#xlims!(1.65, 1.815)
#ylims!(1.73, 1.77)
display(GLMakie.Screen(), f);

#= f = Figure(size = (1000, 600))
ax = Axis(f[1, 1], xgridvisible = false, ygridvisible = false,
xlabel = L"x_1", ylabel = L"x_2", xlabelsize = labelsize, ylabelsize = labelsize,
xticklabelsize = ticksize, yticklabelsize = ticksize)#,
#xticks = [1.786, 1.788, 1.79], yticks = [1.7407, 1.74075, 1.7408], ytickformat = "{:.5f}")
scatter!(tr[tstartpo:tendpo, 1], tr[tstartpo:tendpo, 3], color = :red, markersize = 1.0)
#xlims!(1.7954, 1.7956)
display(GLMakie.Screen(), f); =#