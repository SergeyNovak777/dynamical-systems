username = "sergey"
pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, DifferentialEquations, DynamicalSystems

function rulkov_two_coupled_chem_mem(u, p, t)

  function right_part_x(x, y, z, α)
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

  function right_part_y(x, y, sum_inhibitory, μ, σ)
      return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
  end

  function xi(x, x_th)
      if x > x_th
          return 1.0;
      elseif x <=x_th
          return 0.0;
      end
  end

  ρ(L, k1 ,k2) = k1 + k2 * L^2;

  x1, y1, z1, x2, y2, z2, L = u
  α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2, k1, k2 = p

  if g1 > 0 && g2 >0
      k = 2
  else
      k = 1
  end

  I21 = g2 * ( x_rp - x1 ) * xi(x2, x_th)
  x1n = right_part_x(x1, y1 + (β_syn/k) * I21 + ρ(L, k1, k2) * (x2 - x1), z1, α)
  y1n = right_part_y(x1, y1, (σ_syn/k) * I21, μ, σ)
  z1n = x1

  I12 = g1 * ( x_rp - x2 ) * xi(x1, x_th)
  x2n = right_part_x(x2, y2 + (β_syn/k) * I12 + ρ(L, k1, k2) * (x1 - x2), z2, α)
  y2n = right_part_y(x2, y2, (σ_syn/k) * I12, μ, σ)
  z2n = x2

  L = x1 - x2;
  
  return SVector{7}(x1n, y1n, z1n, x2n, y2n, z2n, L)
end

function get_params_rulkov_two_coupled_chem_mem()
  α = 5.6; σ = 0.2; μ = 0.001;
  β_syn = 0.0001; σ_syn = 1.0; k = 1.0;
  x_rp = -1.5; x_th = -0.8;
  γ_1 = 0.0; γ_2 = 0.0;
  g1 = 0.0; g2 = 0.0;
  k1 = 0.0; k2 = 0.0;
  return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, γ_1, γ_2, g1, g2, k1, k2]
end

tspan = (0, 2_500_000);
t_tr = 1_000_000;
t_window_plot = 1_500_000;

params = get_params_rulkov_two_coupled_chem_mem()

params[1] = 3.9; # α
params[2] = 1.0; # σ
params[10] = 5.0; # g1
params[11] = 1.0; # g2

params[12] = 0.1; # k1
params[13] = 0.00; # k2

u0 = SVector(-2.083638390440308, -3.9302148554862937, -2.0867818075436624,
-1.97137793066347, -3.819163877171352, -1.9745518175123469,
-0.11222999003131551);

ds = DeterministicIteratedMap(rulkov_two_coupled_chem_mem, SVector{7}(u0), params)

x1 = x2 = interval(-3, 3)
y1 = y2 = interval(-2, 2)
z1 = z2 = interval(-3, 3)
L = interval(-3, 3)
box = x1 × x1 × y1 × y2 × z1 × z2 × L

fixedpoints(ds, box, tol = 1e-10, order = 1)