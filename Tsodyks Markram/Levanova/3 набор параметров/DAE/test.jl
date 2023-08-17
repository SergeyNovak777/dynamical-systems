user = "Alex";
pathtorepo = "C:\\Users\\" * user * "\\Desktop\\";

using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\");
include("C:\\Users\\" * user * "\\Desktop\\dynamical-systems\\system.jl")

using StaticArrays, DifferentialEquations, LinearAlgebra, GLMakie;

const τ = 0.013;  const τD = 0.07993;  const τy = 3.3;  const J = 3.07;  const β = 0.300;
const xthr = 0.75; const ythr = 0.4; const α = 1.58; const ΔU0 = 0.305;
M = Diagonal([τ, 1.0, 1.0])

@inbounds function TM_(du, u, p, t)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) )
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) )
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) )
    
    U_ = U(u[3], p)
    du[1] = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) )
    du[2] = (1.0 - u[2]) / p[3] - U_*u[2]*u[1]
    du[3] = (-u[3])/p[4] + p[10] * σ(u[2], p)
    
    nothing
end

function plot3(u, ts, tf)

    indexx,indexy,indexz = 2, 3, 1
    
    GLMakie.activate!()
    f = Figure(resolution = (900, 600))
    axis3 = Axis3(f[1, 1])
    
    lines!(axis3, u[indexx, ts:tf], u[indexy, ts:tf], u[indexz, ts:tf], linewidth = 1.5, color = :black)
    display(GLMakie.Screen(), f)
end

function main()
    tstart = 0.0;
    tfinish = -200.0;
    tspan = (tstart, tfinish);

    I0 = -1.7064; U0 = 0.265;
    p = [α, τ, τD, τy, J, xthr, ythr, U0, ΔU0, β, I0];
    u0 = [8.254587677389917, 0.7099299111721875, 0.549798431121074];

    f = ODEFunction(TM_, mass_matrix = M);
    prob_DAE = ODEProblem(f, u0, tspan, p);
    sol_DAE = solve(prob_DAE, alg = RadauIIA5(), reltol = 1e-8, abstol = 1e-10);

    plot3(sol_DAE, 1, length(sol_DAE));
    return sol_DAE;
end