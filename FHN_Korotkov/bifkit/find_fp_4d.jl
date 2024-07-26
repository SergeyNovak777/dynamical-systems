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
    k1 = 0.09353383; k2 = 64.76190476190476; # 64.76190476190476;
    return (ϵ = ϵ, a = a, g = g, k = k, σ = σ, α = α, k1 = k1, k2 = k2)
end

function FHN2_4d(u, p)
    x1, y1, x2, y2 = u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + (k1 + k2 * (y1 - y2)^2) * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + (k1 + k2 * (y1 - y2)^2) * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    return [dx1dt, dy1dt, dx2dt, dy2dt]
end

u0 = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071]

params = FHN2_try3_params_set();

opt_new = NewtonPar(tol = 1e-9, max_iterations = 10, linsolver=GMRESKrylovKit(), verbose =true);

pmax = 0.15 #0.09353383;
pmin = 0.0 #0.09353383;

opts_con = ContinuationPar(p_min = pmin, p_max = pmax,
                            ds = 0.00001, dsmin = 1e-10, dsmax = 1e-5,
                            nev = 5, detect_bifurcation = 3, newton_options  = opt_new,
                            max_steps  = 1000)


prob =  BifurcationProblem(FHN2_4d, u0, params, (@lens _.k1));

br = continuation(prob, PALC(), opts_con, verbosity=2, linear_algo = BorderingBLS(opt_new.linsolver))

plot(br)


opts_con_k2 = ContinuationPar(p_min = 30.0, p_max = 80.0,
ds = 0.01, dsmin = 1e-10, dsmax = 0.01,
nev = 5, detect_bifurcation = 3, newton_options  = opt_new,
max_steps  = 20000)

hp_codim2_1 = continuation(br, 1, (@lens _.k2),
        opts_con_k2,
        detect_codim2_bifurcation = 2,
        update_minaug_every_step = 1,
        verbosity = 2, 
        linear_algo = BorderingBLS(opt_new.linsolver),
        bothside = true)

plot(hp_codim2_1)