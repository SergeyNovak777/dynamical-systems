user = "Alex";
pathtorepo = "C:\\Users\\" * user * "\\Desktop\\";

using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\");

using StaticArrays, DifferentialEquations, DynamicalSystems, LinearAlgebra, GLMakie;

@inbounds function lorenz(u, p, t)
    x = p[1] * ( u[2] - u[1] );
    y = u[1] * ( p[2] - u[3] ) - u[2];
    z = u[1] * u[2] - p[3] * u[3];
    SVector(x, y, z);
end

function plot3(u, ts, tf)

    indexx,indexy,indexz = 1, 2, 3
    
    GLMakie.activate!()
    f = Figure(resolution = (900, 600))
    axis3 = Axis3(f[1, 1])
    
    lines!(axis3, u[indexx, ts:tf], u[indexy, ts:tf], u[indexz, ts:tf], linewidth = 1.5, color = :black)
    display(GLMakie.Screen(), f)
end

function main()
    σ = 10.0; ρ = 28.0; β = 8/3;
    p = [σ, ρ, β];

    tstart = 0.0; tfinish = 100.0; tspan = (tstart, tfinish);
    u0 = SA[0.001, 0.0, 0.0];


    prob = ODEProblem(lorenz, u0, tspan, p);
    sol = solve(prob, alg = Vern9(), reltol = 1e-8, abstol = 1e-12);

    plot3(sol, 1, length(sol));

end