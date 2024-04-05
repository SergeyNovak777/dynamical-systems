username = "Alex";
env = "integrate";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");

using LinearAlgebra, JLD
using DynamicalSystems, GLMakie, StaticArrays, DifferentialEquations

@inbounds function TM_DAE(du, u, p, t)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) )
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) )
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) )
    
    U_ = U(u[3], p)
    du[1] = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) )
    du[2] = (1.0 - u[2]) / p[3] - U_*u[2]*u[1]
    du[3] = (-u[3])/p[4] + p[10] * σ(u[2], p)
    
    nothing
end

@inbounds function jacob_TM_(u, p, t)
    
    U(y, p, exp50) = p[8] + p[9] / ( 1.0 + exp50 )
    U_y(y, p, exp50) = (50.0 * p[9] * exp50) / (1.0 + exp50)^2
    g(E, x, y, p, U_) = exp((p[5]  * U_ * x * E + p[11]) / p[1])
    σ_der(x, p) = exp( (-20.0) * (x - p[6]) )
    exp50 = exp(-50.0 * (u[3] - p[7]))
    
    U_ = U(u[3], p, exp50)
    Uy = U_y(u[3], p, exp50)
    g_ = g(u[1], u[2], u[3], p, U_)
    σ_deri = σ_der(u[2], p)
    
    g_plus = 1.0 + g_
    g_mult = g_ * U_
    g_plus_mult = p[2] * (g_plus)
    u1p5 = p[5] * u[1]
    Uyu2 = Uy * u[2]
    
    E_E = (-1.0 + ((p[5] * u[2] * g_mult)) / (g_plus) ) / p[2]
    E_x = (u1p5 * g_mult) / (g_plus_mult)
    E_y = (u1p5 * Uyu2 * g_) / (g_plus_mult)
    
    x_E = -U_ * u[2]
    x_x = -1.0 / p[3] - U_ * u[1]
    x_y = -Uyu2 * u[1]
    
    y_x = 20.0 * p[10] * σ_deri / (1.0 + σ_deri)^2
    y_y = -1.0/p[4]
    
    SMatrix{3,3}(E_E, x_E, 0.0,
        E_x, x_x, y_x,
        E_y, x_y, y_y)
end

function get_matrix(fp, p, jac_system, norm = 2)

    Jfp = jac_system(fp, p, 0.0);
    eigen_values_vectors = eigen(Jfp);
    eigen_vectors = eigen_values_vectors.vectors;

    v1 = real(eigen_vectors[:, 1]);
    v2 = real(eigen_vectors[:,2]);
    v3 = imag(eigen_vectors[:, 3]);

    v1 = normalize(v1, norm);
    v2 = normalize(v2, norm);
    v3 = normalize(v3, norm);

    v1 = transpose(v1);
    v2 = transpose(v2);
    v3 = transpose(v3);

    A = transpose([v1; v2; v3]);

    return A;
end

function make_event(fp, ϵ_box, A)
    function condition(u, t, integrator)
        x = Vector(u)
        b = x - fp
        b = Vector(b)
        linprob = LinearProblem(A, b)
        linsolve = solve(linprob)
        return norm(linsolve.u, Inf) - ϵ_box
    end
    return condition
end

affect!(integrator) = terminate!(integrator)

function main()

    p = [1.58, 0.013, 0.07993, 3.3, 3.07, 0.75, 0.4, 0.2650000065028337, 0.305, 0.3, -1.70629651132633]
    M = Diagonal([0.013, 1.0, 1.0])
    fp = [8.345808451620787, 0.7384946646527357, 0.4382985603530577]

    tstart = 0.0; tfinish = 1000.0; tspan = [tstart; tfinish];

    ϵ_box = 1.0e-3;

    A = get_matrix(fp, p, jacob_TM_, 2);
    A = Matrix(A);
    condition = make_event(fp, ϵ_box, A);
    cb = ContinuousCallback(condition, nothing, affect!);

    u0 = [8.344808552431216, 0.7375241488976552, 0.4381235803280678];

    f = ODEFunction(TM_DAE, mass_matrix = M);
    prob_DAE = ODEProblem(f, u0, tspan, p);
    sol_DAE = solve(prob_DAE, GRK4T(), reltol = 1e-8, abstol = 1e-10, callback = cb)

    f = Figure(resolution = (400, 400))
    ax = LScene(f[1, 1], show_axis = true)
    scale!(ax.scene, 50, 50, 1)

    lines!(ax, sol_DAE[2,:], sol_DAE[3, :], sol_DAE[1, :], linewidth = 1.0)
    display(GLMakie.Screen(), f)
    
    return 0;
end