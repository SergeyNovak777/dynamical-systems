username = "irrito"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")


using StaticArrays, DifferentialEquations, DynamicalSystems
using CairoMakie, GLMakie


@inbounds function TM(var, p, t)

    E,x, u = var;
    τ, τ_F, τ_D, U, J, I0, α = p;
    
    g(E, x, u, J, I0) = α * log( 1.0 + exp( (J * u * x * E + I0) / α ))

    du1 = (-E + g(E, x, u, J, I0)) / τ;
    du2 = (1 - x) / τ_D - u * E * x;
    du3 = U * E * (1 - u) - (u - U) / τ_F; 

    return SVector(du1, du2, du3)
end

function TM_params()
    τ = 13.0 / 1000.0
    τ_D = 200.0 / 1000.0
    τ_F = 1500.0 / 1000.0
    J = 3.07 
    U = 0.3
    α = 1.5
    I0 = -1.765
    return [τ, τ_F, τ_D, U, J, I0, α];
end


function get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters)
    integrator_setting = (alg = alg, adaptive = adaptive, abstol = abs_tol, reltol = rel_tol, maxiters = max_iters);
    return integrator_setting;
end

t_truncate(t) = floor(Int64, t / 2)

alg = Vern9()
abs_tol = 1e-10;
rel_tol = 1e-10;
max_iters = 1e8;
adaptive = true;

integrator_setting = get_set_integ_setting(alg, adaptive, abs_tol, rel_tol, max_iters);

parameters = TM_params();

u0_start = [0.208801011940326, 0.47521153582389164, 0.8036140687185943] #[3.44, 0.46, 0.84]
# [7.0, 0.2, 0.7];
# [1.208801011940326, 0.47521153582389164, 0.8036140687185943]
u0_start = SVector{3}(u0_start);

TM(u0_start, parameters, 0)

t_end = 15_000;
tspan = (0.0, t_end);

prob = ODEProblem(TM, u0_start, tspan, parameters)

sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters)

#= ds = CoupledODEs(TM, sol[end], parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 5000);
println("LSE: $(LSE)"); =#

labelsize = 50;
ticksize = 25;

length_sol = length(sol);
ttr = t_truncate(length_sol)
t_plot_start =  ttr
t_plot_end = t_plot_start + 50_000

indexx = 2; indexy  = 3; indexz = 3;
f = Figure(size = (1200 ,600));
ax = Axis3(f[1, 1], xlabel = L"x", ylabel = L"u", zlabel = L"E",
    xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    xgridvisible = false, ygridvisible = false, zgridvisible = false)#,
lines!(ax, sol[indexx, t_plot_start:t_plot_end], sol[indexy, t_plot_start:t_plot_end],
        sol[indexz, t_plot_start:t_plot_end], linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);