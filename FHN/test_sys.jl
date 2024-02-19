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
end

using DifferentialEquations, StaticArrays, GLMakie

# include("/home/sergey/work/repo/dynamical-systems/FHN/test_sys.jl")
function two_FHN(u, p, t)

    x1, y1, x2, y2, z = u
    ϵ, g, k, σ, a, α, k1, k2 = p

    ϕ(y_i, x_i) = atan(y_i / x_i)
    σdiv = σ / 2
    I(ϕ_i) = g / ( 1 + exp( k * ( cos( σdiv ) - cos( ϕ_i - α - σdiv ) ) ) )
    ρz = k1 + k2 * z ^ 2

    dx1dt = ( x1 - x1 ^ 3 / 3 - y1 + I( ϕ(y2, x2) ) + ρz * ( x2 - x1 ) ) / ϵ
    dy1dt = x1 - a
    dzdt = x1 - x2
    dx2dt = ( x2 - x2 ^ 3 / 3 - y2 + I( ϕ(y1, x1) ) + ρz * dzdt ) / ϵ
    dy2dt = x2 - a
    

    return SVector(dx1dt, dy1dt, dx2dt, dy2dt, dzdt)

end

function two_FHN_get_params()
    ϵ = 0.01; a = -1.01
    k = 50.0; g = 0.1; α = 160; σ = 50;
    k1 = 0.09; k2 = 0.0;
    params = [ϵ, g, k, σ, a, α, k1, k2]
    return params
end

function two_FHN_get_helps(param)
    indexparams = "ϵ = 1, g = 2, k = 3, σ = 4, a = 5, α = 6, k1 = 7, k2 = 8";
    keyp = split(nameparam, ", ");
    dict = Dict(zip(keyp, param));

    return dict, indexparams;
end

params = two_FHN_get_params()
u0 = [-1, -0.5, 1, 0.5, 0.6]
tspan = (0.0, 2000.0)

prob = ODEProblem(two_FHN, SVector{length(u0)}(u0), tspan, params)
sol = solve(prob, RK4(), adaptive = false, dt = 0.001, maxiters = 5e6)

tstart = 1; tend = 2000000
f = Figure(size = (1000, 400))
ax = Axis(f[1, 1])
lines!(sol.t[tstart:tend], sol[1, tstart:tend], color = :red, linewidth = 2.0)
lines!(sol.t[tstart:tend], sol[3, tstart:tend], color = :green, linewidth = 2.0)
display(f)

tstart = 1; tend = 2000000
f = Figure(size = (1000, 400))
ax = Axis(f[1, 1])
lines!(sol.t[tstart:tend], sol[2, tstart:tend], color = :red, linewidth = 2.0)
lines!(sol.t[tstart:tend], sol[4, tstart:tend], color = :green, linewidth = 2.0)
display(f)