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
using StaticArrays, DifferentialEquations, DynamicalSystems, CairoMakie, GLMakie, JLD2

function three_coupled_rulkov_first_iteration(u, p)

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

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, x2, y2, z2, x3, y3, z3 = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 = g2 * ( x_rp - x1 ) * xi(x2)
    I31 = g1 * ( x_rp - x1 ) * xi(x3)

    I12 = g1 * ( x_rp - x2 ) * xi(x1)
    I32 = g2 * ( x_rp - x2 ) * xi(x3)

    I13 = g2 * ( x_rp - x3 ) * xi(x1)
    I23 = g1 * ( x_rp - x3 ) * xi(x2)

    x1n = right_part_x(x1, y1 + (β_syn/k) * (I21 + I31), z1)    
    y1n = right_part_y(x1, y1, (σ_syn/k) * (I21 + I31))
    z1n = x1

    x2n = right_part_x(x2, y2 + (β_syn/k) * (I12 + I32), z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * (I12 + I32))
    z2n = x2

    x3n =  right_part_x(x3, y3 + (β_syn/k) * (I13 + I23), z3)
    y3n = right_part_y(x3, y3, (σ_syn/k) * (I13 + I23))
    z3n = x3

    return SVector(x1n, y1n, z1n, I21, I31, x2n, y2n, z2n, I12, I32, x3n, y3n, z3n, I13, I23)
end

function three_coupled_rulkov(u, p, t)

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

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, I21prev, I31prev, x2, y2, z2, I12prev, I32prev, x3, y3, z3, I13prev, I23prev = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2 = p

    if g1 > 0 && g2 >0
        k = 2
    else
        k = 1
    end

    I21 =  γ_2 * I21prev + g2 * ( x_rp - x1 ) * xi(x2)
    I31 =  γ_1 * I31prev + g1 * ( x_rp - x1 ) * xi(x3)
    x1n = right_part_x(x1, y1 + (β_syn/k) * (I21 + I31), z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * (I21 + I31))
    z1n = x1

    I12 = γ_1 * I12prev + g1 * ( x_rp - x2 ) * xi(x1)
    I32 = γ_2 * I32prev + g2 * ( x_rp - x2 ) * xi(x3)
    x2n = right_part_x(x2, y2 + (β_syn/k) * (I12 + I32), z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * (I12 + I32))
    z2n = x2

    I13 = γ_2 * I13prev +  g2 * ( x_rp - x3 ) * xi(x1)
    I23 = γ_1 * I23prev + g1 * ( x_rp - x3 ) * xi(x2)
    x3n =  right_part_x(x3, y3 + (β_syn/k) * (I13 + I23), z3)
    y3n = right_part_y(x3, y3, (σ_syn/k) * (I13 + I23))
    z3n = x3

    return SVector(x1n, y1n, z1n, I21, I31, x2n, y2n, z2n, I12, I32, x3n, y3n, z3n, I13, I23)
end

function get_params_three_coupled_rulkov()
    α = 3.9; σ =1.0; μ = 0.001;
    β_syn = 0.0001; σ_syn = 1.0;
    x_rp = -1.5; x_th = -0.8;
    γ_1 = 0.0; γ_2 = 0.0;
    g1 = 1.0; g2 = 7.0;
    return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2]
end

Hs(x, k) = Statistics.mean(x) + k * Statistics.std(x)

function CALCPDF_debil(spikes, threshold, ϵ)
    ee_counter = [sum(i->s<=i<s+ϵ, spikes) for s in threshold]
    pdf = ee_counter ./ length(spikes)
    return pdf
end

params = get_params_three_coupled_rulkov()
#params[8] = 0.5;
#params[9] = 0.5;
params[10] = 0.0;
params[11] = 1.0;
tspan = (0, 100000);

u0 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
u0_first_iteration = three_coupled_rulkov_first_iteration(u0, params);
prob = DiscreteProblem(three_coupled_rulkov, SVector{length(u0_first_iteration)}(u0_first_iteration), tspan, params);
sol = solve(prob);
Ttr = 50_000;
point_from_attractor = sol[:, Ttr]

ds = DeterministicIteratedMap(three_coupled_rulkov, point_from_attractor, params)

index_control_parameter = 10;
index_save_variable = 1 #[1,6,11];
g1_range = range(0.0, 1.0, length = 2000);
n = 1000;
Ttr = 5000;

output = orbitdiagram(ds, index_save_variable, index_control_parameter, g1_range; n, Ttr)

L = length(g1_range)
x = zeros(n*L)
x1 = zeros(n*L)
#= x2 = copy(x)
x3 = copy(x) =#
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= g1_range[j]
    #first_elements = [v[1] for v in output[j]]
    #second_elements = [v[2] for v in output[j]]
    #third_elements = [v[3] for v in output[j]]
    x1[(1 + (j-1)*n):j*n] .=  output[j] # first_elements
    #x2[(1 + (j-1)*n):j*n] .= second_elements #output[j]
    #x3[(1 + (j-1)*n):j*n] .= third_elements #output[j]
end

f = Figure();
ax = Axis(f[1, 1], ylabel = L"x_1")
scatter!(ax, x, x1, markersize = 0.5, color = :black)
display(GLMakie.Screen(), f)

#= f = Figure();
ax = Axis(f[1, 1], ylabel = L"x_2")
scatter!(ax, x, x2, markersize = 0.5, color = :black)
display(GLMakie.Screen(), f)

f = Figure();
ax = Axis(f[1, 1], ylabel = L"x_3")
scatter!(ax, x, x3, markersize = 0.5, color = :black)
display(GLMakie.Screen(), f)

xsum = x1 + x2 + x3

f = Figure();
ax = Axis(f[1, 1], ylabel = L"x_{sum}")
scatter!(ax, x, xsum, markersize = 0.5, color = :black)
display(GLMakie.Screen(), f) =#