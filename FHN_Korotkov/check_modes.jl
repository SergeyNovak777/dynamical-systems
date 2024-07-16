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

function conv_sol_u(sol_u)
    len_sol = length(sol_u);
    dim = length(sol_u[1]);
    sol = zeros(len_sol, dim);

    for index in range(1, dim)
        sol[:, index] = [v[index] for v in sol_u];
    end

    return sol;
end

function get_sol(prob, integrator_setting)

    sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol,
                maxiters = integrator_setting.maxiters);
    sol_t = sol.t;
    sol_u = conv_sol_u(sol.u);
    return sol_t, sol_u;
end

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9();
adaptive = true;
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = FHN2_try3_params();
parameters[7] = 0.096859
parameters[8] = 4.448160535117057

u0_start = [-1.002201845684165, -0.6387935561024373, -1.0254528312514828, -0.6308174333307994, -0.007976122771661428]
u0_start = SVector{5}(u0_start);

t_end = 20_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_try3, u0_start, tspan, parameters)
sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol,
                maxiters = integrator_setting.maxiters);

ds = CoupledODEs(FHN2_try3, u0_start, parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 10000);
println(LSE);

len_sol = length(sol.t);
ttr = t_truncate(len_sol);

labelsize = 40;
ticksize = 30;
t_plot_start = ttr;
t_plot_end = len_sol;

f = Figure();
ax = Axis3(f[1, 1]);
lines!(ax, sol[1, t_plot_start:t_plot_end], sol[3, t_plot_start:t_plot_end],
        sol[2, t_plot_start:t_plot_end], linewidth = 0.5, color = :black);
display(GLMakie.Screen(), f)