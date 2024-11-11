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
parameters[3] = 0.01;
parameters[7] = 0.0; #0.09;
parameters[8] = 0;   ; # 75.7;

u0_start = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];
u0_start = SVector{4}(u0_start);

t_end = 5_000;
tspan = (0.0, t_end);

prob = ODEProblem(FHN2_4d, u0_start, tspan, parameters)

sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
                abstol = integrator_setting.abstol, reltol = integrator_setting.reltol, 
                maxiters = integrator_setting.maxiters);


#= function FHN2_try3_params()
    ϵ = 0.01; a = -1.01;
    g = 0.1; k = 50.0; σ = 50.0 * pi / 180; α = 160.0 * pi / 180;
    k1 = 0.0; k2 = 0.0
    return [ ϵ, a, g, k, σ, α, k1, k2]
end =#


Ttr = 150_000;

I(ϕ_i, g, k, σ, α) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))

array_ϕ_i = atan.(sol[2, Ttr:end], sol[1, Ttr:end])
array_I = I.(array_ϕ_i, parameters[3], parameters[4], parameters[5], parameters[6])

labelsize = 40;
ticksize = 25       ;
CairoMakie.activate!();

                
f = Figure(size = (400 ,400));
ax = Axis(f[1, 1], xlabel = L"ϕ", ylabel = L"I", 
    xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize)

lines!(ax,array_ϕ_i, array_I,
        linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);

f = Figure(size = (400 ,400));
ax = Axis(f[1, 1], xlabel = L"x_1", ylabel = L"y_1", 
    xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize)

lines!(ax, sol[1, Ttr:end], sol[2, Ttr:end], linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);


t_series_t_plt_end = Ttr + 15_000;
f = Figure(size = (400 ,400));
ax = Axis(f[1, 1], xlabel = L"t", ylabel = L"x_1", 
    xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = ticksize, yticklabelsize = ticksize)

lines!(ax, sol.t[Ttr:t_series_t_plt_end], sol[1, Ttr:t_series_t_plt_end], linewidth = 1.5, color = :black);
display(GLMakie.Screen(), f);   


ds = CoupledODEs(FHN2_4d, sol[end], parameters,
diffeq = integrator_setting);

LSE = lyapunovspectrum(ds, 10_000);
println("LSE: $(LSE)");