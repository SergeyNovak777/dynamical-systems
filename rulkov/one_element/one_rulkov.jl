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

using StaticArrays, Statistics, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie

function two_coupled_rulkov(u, p, t)

    function right_part_x(x, y, z)
        if x <= 0.0
            return α/(1.0 - x) + y
        elseif (0.0 < x < (α + y)) && (z <=0.0)
            return α  +y
        elseif (x >= α + y) || (z > 0)
            return -1.0
        else
            return -1.0;
        end
    end

    function right_part_y(x, y)
        return y + μ * ( -x - 1.0 + σ )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x, y, z = u
    α, σ, μ = p

    xn = right_part_x(x, y, z)
    yn = right_part_y(x, y)
    zn = x

    return SVector{3}(xn, yn, zn)
end


function get_params_two_coupled_rulkov()
    α = 4.6; σ =-0.1; μ = 0.001;
    return [α, σ, μ]
end

params = get_params_two_coupled_rulkov();
params[1] = 3.5;#4.4;
params[2] = 0.15;#0.1;
u0 = [1.0, 0.3, 0.01];
tspan = (0, 100000);

prob = DiscreteProblem(two_coupled_rulkov, SVector{3}(u0), tspan, params);
sol = solve(prob);

t_start = 50000;
t_end = 80000;

f = Figure()
ax = Axis(f[1, 1])
lines!(ax, sol.t[t_start:t_end], sol[1, t_start:t_end], linewidth = 1.0, color = :blue)
display(GLMakie.Screen(), f)

f = Figure()
ax = Axis(f[1, 1])
scatter!(ax, sol[1, t_start:t_end], sol[2, t_start:t_end], markersize = 1.0, color = :blue)
display(GLMakie.Screen(), f)

ds = DeterministicIteratedMap(two_coupled_rulkov, sol[end], params)
Λs = lyapunovspectrum(ds, 500000)

