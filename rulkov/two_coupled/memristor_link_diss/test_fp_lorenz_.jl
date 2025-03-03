username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie
using BenchmarkTools
function lorenz_rule(u, p, t)
    σ = p[1]; ρ = p[2]; β = p[3]
    du1 = σ*(u[2]-u[1])
    du2 = u[1]*(ρ-u[3]) - u[2]
    du3 = u[1]*u[2] - β*u[3]
    return SVector{3}(du1, du2, du3)
end
function lorenz_jacob(u, p, t)
    σ, ρ, β = p
    return SMatrix{3,3}(-σ, ρ - u[3], u[2], σ, -1, u[1], 0, -u[1], -β)
end

ρ, β = 30.0, 10/3
lorenz = CoupledODEs(lorenz_rule, 10ones(3), [10.0, ρ, β])
# Define the box within which to find fixed points:
x = y = interval(-20, 20)
z = interval(0, 40)
box = x × y × z

#fp, eigs, stable = 


@time fixedpoints(lorenz, box)

@time fixedpoints(lorenz, box, lorenz_jacob)