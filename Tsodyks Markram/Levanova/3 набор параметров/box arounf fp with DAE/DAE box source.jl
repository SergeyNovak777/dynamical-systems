username = "Alex";
env = "integrate";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");

include(pathtorepo * "dynamical-systems\\system.jl");
include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\box arounf fp with DAE\\DAE box header.jl")

using LinearAlgebra, JLD
using DynamicalSystems, GLMakie, StaticArrays, DifferentialEquations

function main(number_point)

    # load fixed parameters and mass matrix
    p, M = load_params()

    # load homoclinic curve
    I0_hom, U0_hom = load_hom_curve()

    # choose value for control parameters
    param_index_from_hom_curve = 1;
    I0 = I0_hom[param_index_from_hom_curve];
    U0 = U0_hom[param_index_from_hom_curve];
    p[8] = U0; p[11] = I0;
    
    # time of intregrate
    tstart = 0.0; tfinish = 1000.0; tspan = [tstart; tfinish];

    # size of box
    ϵ_box = 1.0e-5;

    # interval for search fixed points
    E = interval(-40, 40); x = interval(-5, 5); y = x;
    # space for search fixed points
    box = IntervalBox(E, x, y);

    # counts points on side box
    count_points_on_side = 21;

    # find fixed points
    fp = get_fp(TM, jacob_TM_, p, box, 3);
    # check fixed points
    fp = check_fp(fp);
    fp = Vector(fp);

    # get matrix from eigenvectors
    A = get_matrix(fp, p, jacob_TM_, 2)
    A = Matrix(A);
    left_side, right_side = get_left_right_side(ϵ_box, count_points_on_side, fp, A, 3);
    up_side, down_side = get_up_down_side(ϵ_box, count_points_on_side, fp, A, 3);

    points = cat(left_side, right_side, up_side, down_side,  dims = 1);

    condition = make_event(fp, ϵ_box, A);
    cb = ContinuousCallback(condition, nothing, affect!);
    #number_point = 7;
    u0 = [left_side[number_point, 1], left_side[number_point, 2], left_side[number_point, 3]];

    f = ODEFunction(TM_DAE, mass_matrix = M);
    prob_DAE = ODEProblem(f, u0, tspan, p);
    sol_DAE = solve(prob_DAE, GRK4T(), reltol = 1e-10, abstol = 1e-14, callback = cb);

    f = Figure(resolution = (400, 400))
    ax = LScene(f[1, 1], show_axis = true)
    scale!(ax.scene, 50, 50, 1)

    lines!(ax, sol_DAE[2,:], sol_DAE[3, :], sol_DAE[1, :], linewidth = 1.0);
    scatter!(ax, sol_DAE[2,1], sol_DAE[3, 1], sol_DAE[1, 1], markersize = 5.0, color = :green);
    scatter!(ax, sol_DAE[2,end], sol_DAE[3, end], sol_DAE[1, end], markersize = 5.0, color = :red);

    scatter!(ax, left_side[:, 2], left_side[:, 3], left_side[:, 1], markersize = 5.0, color = :black);
    scatter!(ax, right_side[:, 2], right_side[:, 3], right_side[:, 1], markersize = 5.0, color = :black);
    scatter!(ax, up_side[:, 2], up_side[:, 3], up_side[:, 1], markersize = 5.0, color = :blue);
    scatter!(ax, down_side[:, 2], down_side[:, 3], down_side[:, 1], markersize = 5.0, color = :blue);
    display(GLMakie.Screen(), f)

    return sol_DAE;
end