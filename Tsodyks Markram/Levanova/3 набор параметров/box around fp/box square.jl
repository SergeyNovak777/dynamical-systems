const τ = 0.013;  const τD = 0.07993;  const τy = 3.3;  const J = 3.07;  const β = 0.300
const xthr = 0.75; const ythr = 0.4; const α = 1.58; const ΔU0 = 0.305

username = "Alex";
env = "integrate";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");
include(pathtorepo * "dynamical-systems\\system.jl");
include(pathtorepo * "dynamical-systems\\header.jl");

using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2,
GLMakie, LinearAlgebra, LinearSolve;

homoclinic_curve = load_hom_curve();
index_point_from_curve = 1;

I0 = homoclinic_curve[1][index_point_from_curve];
U0 = homoclinic_curve[2][index_point_from_curve];
p = SA[α, τ, τD, τy, J, xthr, ythr, U0, ΔU0, β, I0];

tstart = 0.0; tfinish = 1000.0; tspan = [tstart; tfinish];
ϵ_box = 1.0e-5; ϵ_shift = 1.0e-4;
E = interval(-40, 40); x = interval(-5, 5); y = x;
box = IntervalBox(E, x, y);

count_points_on_side = 11;

fp = get_fp(TM, jacob_TM_, p, box, 3);
if typeof(fp) == StateSpaceSet{3, Float64}
    fp = fp[1]
else
    println("more one fixed point")
end

A = get_matrix(fp, p, jacob_TM_, Inf);

points_left_side, points_right_side = left_right_side(ϵ_box, count_points_on_side, fp, A, 3);
points_up_side, points_down_side = up_down_side(ϵ_box, count_points_on_side, fp, A, 3);

points = cat(points_left_side, points_right_side, points_up_side, points_down_side,  dims = 1);

#norms, αs = get_norm_αs_(points, fp, A);

condition = make_event(fp, ϵ_box, A);
cb = ContinuousCallback(condition, nothing, affect!);
u0 = SA[points_left_side[1, 1], points_left_side[1, 2], points_left_side[1, 3]];

prob = ODEProblem(TM, u0, tspan, p);
sol = solve(prob, alg = Vern9(), abstol = 1e-12, reltol = 1e-10, callback = cb);

