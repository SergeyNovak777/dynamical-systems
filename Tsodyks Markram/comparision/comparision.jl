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

using DifferentialEquations, StaticArrays, CairoMakie, GLMakie, StaticArrays
t_truncate(t) = floor(Int64, t / 2)


function TM_4_get_params_comp()

    τ = 0.013; τ_D = 0.07993; τ_y = 3.3; τ_F = 0.03;
    α = 1.58; β = 0.300; J = 3.07;
    ΔU0 = 0.305;
    I0 = -1.8; U0 = 0.265
    xthr = 0.75; ythr = 0.4

    return [τ, τ_D, τ_F, τ_y, α, β, J, ΔU0, xthr, ythr, U0, I0]
end

function TM_model_get_params_comp()

    τ = 0.013; τD = 0.07993; τy = 3.3; J = 3.07; β = 0.300
    xthr = 0.75; ythr = 0.4
    α = 1.58; ΔU0 = 0.305;
    # control parameters
    I0 = -1.8; U0 = 0.265
    params = [α, τ, τD, τy, J, xthr, ythr, U0, ΔU0, β, I0]

    return params
end


@inbounds function TM(u, p, t)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) )
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) )
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) )
    
    U_ = U(u[3], p)
    du1 = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) ) / p[2]
    du2 = (1.0 - u[2]) / p[3] - U_*u[2]*u[1]
    du3 = (-u[3])/p[4] + p[10] * σ(u[2], p)
    
    return SVector(du1, du2, du3)
end

α, τ, τD, τy, J, xthr, ythr, U0, ΔU0, β, I0 = TM_model_get_params_comp()

sys_TM3 = TM;
params_TM3 = TM_model_get_params_comp();
u0_TM3 = [0.6677760562215335, 0.9635510440831664, 0.4484567117090125]

U(y) = U0 + ΔU0 / ( 1 + exp( -50.0 * ( y - ythr ) ) )

sys_TM4 = TM_4;
params_TM4 = TM_4_get_params_comp();

u0_TM4 = [0.6677760562215335, 0.9635510440831664, U(u0_TM3[3]) , 0.4484567117090125]

tspan = (0.0, 2_000.0);

prob_TM4 = ODEProblem(sys_TM4, u0_TM4, tspan, params_TM4)
sol_TM4 = solve(prob_TM4, RK4(),  adaptive = false, dt = 0.01, maxiters = 5e6);


prob_TM3 = ODEProblem(sys_TM3, u0_TM3, tspan, params_TM3);
sol_TM3 = solve(prob_TM3, RK4(),  adaptive = false, dt = 0.01, maxiters = 5e6);


length_sol = length(sol_TM4);
Ttr = t_truncate(length_sol);


t_plot_start = 1; t_plot_end = 10_000;
width_windown = 1600;
hight_window = 600;
lw = 2.0;

f = Figure(size = (width_windown, hight_window))
ax_E = Axis(f[1, 1], xlabel = "t", ylabel = "E", xlabelsize = 22, ylabelsize = 22,xticklabelsize = 17, yticklabelsize = 17)
ax_X = Axis(f[2, 1], xlabel = "t", ylabel = "X", xlabelsize = 22, ylabelsize = 22,xticklabelsize = 17, yticklabelsize = 17)
ax_u = Axis(f[3, 1], xlabel = "t", ylabel = "u", xlabelsize = 22, ylabelsize = 22,xticklabelsize = 17, yticklabelsize = 17)
ax_y = Axis(f[4, 1], xlabel = "t", ylabel = "Y", xlabelsize = 22, ylabelsize = 22,xticklabelsize = 17, yticklabelsize = 17)

lines!(ax_E, sol_TM4.t[t_plot_start:t_plot_end], sol_TM4[1,t_plot_start:t_plot_end], color = "blue", linewidth = lw)
lines!(ax_X, sol_TM4.t[t_plot_start:t_plot_end], sol_TM4[2,t_plot_start:t_plot_end], color = "blue", linewidth = lw)
lines!(ax_u, sol_TM4.t[t_plot_start:t_plot_end], sol_TM4[3,t_plot_start:t_plot_end], color = "blue", linewidth = lw)
lines!(ax_y, sol_TM4.t[t_plot_start:t_plot_end], sol_TM4[4,t_plot_start:t_plot_end], color = "blue", linewidth = lw)


lines!(ax_E, sol_TM3.t[t_plot_start:t_plot_end], sol_TM3[1,t_plot_start:t_plot_end], color = "green", linewidth = lw)
lines!(ax_X, sol_TM3.t[t_plot_start:t_plot_end], sol_TM3[2,t_plot_start:t_plot_end], color = "green", linewidth = lw)
lines!(ax_y, sol_TM3.t[t_plot_start:t_plot_end], sol_TM3[3,t_plot_start:t_plot_end], color = "green", linewidth = lw)
display(GLMakie.Screen(), f)


f = Figure(size = (600, 600))
ax = Axis3(f[1, 1], xlabel = L"x", ylabel = L"u", zlabel = L"E")
lines!(sol_TM4[2,t_plot_start:t_plot_end], sol_TM4[4,t_plot_start:t_plot_end], sol_TM4[1,t_plot_start:t_plot_end], color = "blue")
lines!(sol_TM3[2,t_plot_start:t_plot_end], sol_TM3[3,t_plot_start:t_plot_end], sol_TM3[1,t_plot_start:t_plot_end],color = "green")
display(GLMakie.Screen(), f)


