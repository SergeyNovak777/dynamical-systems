pathtorepo = "/home/irrito/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/bifurcation/")

using Revise, Plots;
using BifurcationKit;

function FHN2_4d!(du, u, p, t = 0)
    x1, y1, x2, y2 = u
    (;ϵ, a, g, k, σ, α, k1, k2) = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    du[1] = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + (k1 + k2 * (y1 - y2)^2) * (x2 - x1) ) / ϵ
    du[2] = x1 - a
    du[3] = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + (k1 + k2 * (y1 - y2)^2) * (x1 - x2) ) / ϵ
    du[4] = x2 - a
    
    du
end

#= function FHN2_try3_params_bk_version()
    ϵ = 0.01; a = -1.01;
    g = 0.1; k = 50.0; σ = 50.0 * pi / 180; α = 160.0 * pi / 180;
    k1 = 0.0; k2 = 0.0
    return ( ϵ, a, g, k, σ, α, k1, k2)
end =#

u0 = [-0.9816946043747945, -0.6320919525134647, -1.0342265829731392, -0.638226338524071];

params = (ϵ = 0.01, a = -1.01, g = 0.1, k = 50.0,
            σ = 50.0 * pi / 180, α = 60.0 * pi / 180, k1 = 0.0, k2 = 40.0);

recordFromSolution(x, p; k...) = (x1 = x[1], y1 = x[2], x2 = x[3], y2 = x[4]);

prob = BifurcationProblem(FHN2_4d!, u0, params,
	(@optic _.k1), record_from_solution = recordFromSolution);

opts_br = ContinuationPar(p_min = 0.0, p_max = 0.1, dsmax = 0.01, detect_bifurcation  = 3 , max_steps = 2000);

br = continuation(prob, PALC(), opts_br)

scene = plot(br, legend=:topleft)