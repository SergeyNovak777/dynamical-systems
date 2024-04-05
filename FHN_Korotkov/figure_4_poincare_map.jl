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

# INCLUDE
#---------------------------------------------------------------------
include("/home/sergey/work/repo/dynamical-systems/system.jl")
#----------------------------------------------------------------------
#PACKAGES
using StaticArrays, DifferentialEquations, DynamicalSystems
using CairoMakie
using JLD2
#----------------------------------------------------------------------

function plot_timeseries(sol, tstart, tend, index, labels)
    f = Figure(size = (1000, 250))
    ax = Axis(f[1,1], xlabel = "time", ylabel = labels)
    lines!(ax, sol.t[tstart:tend], sol[index, tstart:tend], color = :black)
    display(f)
end

function plot_two_timeseries(sol, tstart, tend, indexs, labels)
    indexx, indexy = indexs
    f = Figure(size = (1000, 250))
    ax = Axis(f[1,1], xlabel = "time", ylabel = labels)
    lines!(ax, sol.t[tstart:tend], sol[indexx, tstart:tend], color = :red)
    lines!(ax, sol.t[tstart:tend], sol[indexy, tstart:tend], color = :green)
    display(f)
end

function get_u0(x1, y1, x2, y2)
    z = y1 - y2
    return SVector{5}([x1, y1, x2, y2, z])
end

function main(x1, y1, x2, y2)
    parameters = FHN2_try3_params()
    parameters[7] = 0.09
    parameters[8] = 0.0

    tspan = (0.0, 1000)

    integrator_setting = (alg = DP8(), adaptive = true, abstol = 1e-11, reltol = 1e-11)

    u0 = get_u0(x1, y1, x2, y2) # 1.0,2.0,3.0,4.0
    prob = ODEProblem(FHN2_try3, u0, tspan, parameters)
    sol = solve(prob, integrator_setting.alg, adaptive = integrator_setting.adaptive,
    abstol = integrator_setting.abstol, reltol = integrator_setting.abstol)

    len_sol = length(sol.u[:,1])
    tstart = floor(Int64, len_sol / 2); tend = len_sol

    indexs_two_ts = [1, 3] # x1, x2
    label_two_ts = "x_i"
    plot_two_timeseries(sol, tstart, tend, indexs_two_ts, label_two_ts)

    indexs_ts = 5
    label_ts = "z"
    plot_timeseries(sol, tstart, tend, indexs_ts, label_ts)

    integ_set = (alg = Vern9(), adaptive = true, abstol = 1e-13, reltol = 1e-13)
    ds = CoupledODEs(FHN2_try3, sol[end], parameters, diffeq = integ_set)
    pmap = PoincareMap(ds, (4, 0.0))

    tr, trange = trajectory(pmap, 1000000)

    tstartpo = 10000; tendpo= 1000000
    
    f = Figure(size = (450, 450))
    ax = Axis(f[1, 1],
    xlabel = L"x_1", ylabel = L"x_2")
    scatter!(tr[tstartpo:tendpo, 1], tr[tstartpo:tendpo, 3], color = :red, markersize = 1.0)
    xlims!(1.645, 1.85)
    ylims!(1.725, 1.765)
    display(f)


    f = Figure(size = (450, 450))
    ax = Axis(f[1, 1], xgridvisible = false, ygridvisible = false,
    xlabel = L"x_1", ylabel = L"x_2")

    scatter!(tr[tstartpo:tendpo, 1], tr[tstartpo:tendpo, 3], color = :red, markersize = 1.0)
    xlims!(1.7954, 1.7956)
    ylims!(1.741116, 1.741128)
    display(f)

    #jldsave("test_1e6.jld2"; x1= tr[:,1], x2 = tr[:,3])
end

#=
az = -0.3pi
f = Figure(size = (900, 600))
ax3d = Axis3(f[1, 1], xlabel = L"x_1", ylabel = L"x_2", zlabel = L"y_2", xlabelsize = labelsize, ylabelsize = labelsize, zlabelsize = labelsize,
    xticks = [-2, 0, 2], yticks = [-2, 0, 2], xticklabelsize = ticksize, yticklabelsize = ticksize, zticklabelsize = ticksize,
    azimuth = az, xgridvisible = false, ygridvisible = false, zgridvisible = false)
lines!(ax3d, sol[1, tstart:tend], sol[3, tstart:tend], sol[indexz, tstart:tend], linewidth = lw)
display(f)
=#

function test_pmap()
    function duffing_rule(u,p,t)
        d, a, ω = p
        du1 =  u[2]
        du2 =  -u[1] - u[1]*u[1]*u[1] - d*u[2] + a*sin(ω*t)
        return SVector(du1, du2)
    end
    T0 = 25.0
    p0 = [0.1, 7, 2π/T0]
    u0 = [1.1, 1.1]
    ds = CoupledODEs(duffing_rule, u0, p0)
    duffing = StroboscopicMap(ds, T0)
    
    # We want to change both the parameter `ω`, but also the
    # period of the stroboscopic map. `orbitdiagram` allows this!
    Trange = range(8, 26; length = 201)
    ωrange = @. 2π / Trange
    n = 200
    output = orbitdiagram(duffing, 1, 3, ωrange; n, u0, Ttr = 100, periods = Trange)
    
    L = length(Trange)
    x = Vector{Float64}(undef, n*L)
    y = copy(x)
    for j in 1:L
        x[(1 + (j-1)*n):j*n] .= Trange[j]
        y[(1 + (j-1)*n):j*n] .= output[j]
    end
    
    fig, ax = scatter(x, y; axis = (xlabel = L"T", ylabel = L"u_1"),
        markersize = 8, color = ("blue", 0.25),
    )
    ylims!(ax, -1, 1)
    display(fig)
    return ds
end