if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/bifurcation/")
end

using BifurcationKit, Setfield, LinearAlgebra, Plots, Parameters

function FHN2_try3_params_set()
    ϵ = 0.01; a = -1.01;
    g = 0.1; k = 50.0; σ = 50.0 * pi / 180; α = 160.0 * pi / 180;
    k1 = 0.0; k2 = 40.0;
    return (ϵ = ϵ, a = a, g = g, k = k, σ = σ, α = α, k1 = k1, k2 = k2)
end

function FHN2_try3(u, p)
    x1, y1, x2, y2, z= u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))
    ρz = k1 + k2 * z ^ 2

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + ρz * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + ρz * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    dzdt = x1 - x2
    return [dx1dt, dy1dt, dx2dt, dy2dt, dzdt]
end


u0 = [-1.01
        -0.6544965895143954
        -1.01
        -0.6544965895143952
        0.0]#-1.3089931790287905];

params = FHN2_try3_params_set();

opt_new = NewtonPar(tol = 1e-9, max_iterations = 10, linsolver=GMRESKrylovKit(), verbose =true);

pmax = 0.1;
pmin = 0.0;

opts_con = ContinuationPar(p_min = pmin, p_max = pmax,
                            ds = 0.0001, dsmin = 1e-5, dsmax = 0.001,
                            nev = 3, detect_bifurcation = 3, newton_options  = opt_new,
                            max_steps  = 1000)


prob =  BifurcationProblem(FHN2_try3, u0, params, (@lens _.k1));

br = continuation(prob, PALC(), opts_con, verbosity=2, linear_algo = BorderingBLS(opt_new.linsolver))

plot(br)