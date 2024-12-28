username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, DifferentialEquations, DynamicalSystems
using CairoMakie

function standardmap_rule(x, k, n)
    theta = x[1]; p = x[2]
    p += k[1]*sin(theta)
    theta += p
    return SVector(mod2pi(theta), mod2pi(p))
end

standardmap = DeterministicIteratedMap(standardmap_rule, rand(2), [1.0])
xs = range(0, stop = 2π, length = 11); ys = copy(xs)
ics = [SVector{2}(x,y) for x in xs for y in ys]

# All permutations of [±1, ±1]:
singss = lambdaperms(2)[2] # second entry are the signs

# I know from personal research I only need this `inds`:
indss = [[1,2]] # <- must be container of vectors!

λs = 0.005 # <- only this allowed to not be vector (could also be vector)

orders = [2, 3, 4, 5, 6, 8]

ALLFP = []

for o in orders
    FP = periodicorbits(standardmap, o, ics, λs, indss, singss)
    push!(ALLFP, FP)
end

iters = 1000
dataset = trajectory(standardmap, iters)[1]
for x in xs
    for y in ys
        append!(dataset, trajectory(standardmap, iters, [x, y])[1])
    end
end

fig = Figure()
ax = Axis(fig[1,1]; xlabel = L"\theta", ylabel = L"p",
    limits = ((xs[1],xs[end]), (xs[1],xs[end]))
)
scatter!(ax, dataset[:, 1], dataset[:, 2]; markersize = 1, color = "black")
fig


markers = [:diamond, :utriangle, :rect, :pentagon, :hexagon, :circle]

for i in 1:6
    FP = ALLFP[i]
    o = orders[i]
    scatter!(ax, columns(FP)...; marker=markers[i], color = Cycled(i),
        markersize = 30 - 2i, strokecolor = "grey", strokewidth = 1, label = "order $o"
    )
end
axislegend(ax)
fig