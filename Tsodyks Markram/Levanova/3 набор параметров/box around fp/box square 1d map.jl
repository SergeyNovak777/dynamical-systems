const τ = 0.013;  const τD = 0.07993;  const τy = 3.3;  const J = 3.07;  const β = 0.300
const xthr = 0.75; const ythr = 0.4; const α = 1.58; const ΔU0 = 0.305

username = "Alex";
env = "integrate";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
include(pathtorepo * "dynamical-systems\\system.jl");
include(pathtorepo * "dynamical-systems\\header.jl");

using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");
using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2, GLMakie, LinearAlgebra, LinearSolve;

function square(sys, jac_sys, p, tspan,  ϵ_box, box,  counts_number_on_side)
    
    fp = get_fp(sys, jac_sys, p, box, 3);
    if typeof(fp) == StateSpaceSet{3, Float64}
        fp = fp[1];
    else
        println("more one fixed point");
        exit();
    end

    A = get_matrix(fp, p, jacob_TM_, Inf);
    points_up_side, points_down_side = up_down_side(ϵ_box, counts_number_on_side, fp, A, 3);
    points = cat(points_up_side, points_down_side,  dims = 1);

    condition = make_event(fp, ϵ_box, A);
    cb = ContinuousCallback(condition, nothing, affect!);
    println(points);

    """for index in range(1, length(points), step = 1)
        u0 = SA[points[index, 1], points[index, 2], points[index, 3]]
        prob = ODEProblem(sys, u0, tspan, p)

    end"""




end

function main()
    
    # setting for range param
    ϵ_shift = 1.0e-3;
    index_point_from_hom = 1;
    p_len = 5;
    direction = "increase";
    pcontrolname = "I0";

    # time for integrator
    tstart = 0.0;
    tfinish = 1000.0;
    tspan = (tstart, tfinish);

    # setting for search fixed points
    E= interval(-40, 40);
    x = interval(-5, 5);
    y = x;
    box = IntervalBox(E, x, y);

    # setting for square
    ϵ_box = 1e-5;
    counts_number_on_side = 11;
    sys = TM;
    jac_sys = jacob_TM_;

    p_range, pfix = get_range(ϵ_shift, index_point_from_hom, p_len; param = pcontrolname, direction_shift = direction)
   
    for index in range(1, p_len, step = 1)
        
        p = get_set_p(p_range[index], pfix, pcontrolname);
        println(p);
        square(sys, jac_sys, p, tspan,  ϵ_box, box,  counts_number_on_side)
    end

end