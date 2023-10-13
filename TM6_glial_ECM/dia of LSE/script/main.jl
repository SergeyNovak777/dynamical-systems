username = "Alex";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\");

using StaticArrays, DifferentialEquations, DynamicalSystems;

include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\system.jl");


#=
arguments of diagram_lse
* ds
* parameters
* initial condition
* setting of integrator
* index control parameter
* range of control parameter
* name of saving files
=#

function diagram_lse(ds, params, u0, integ_set,
    index_cont_param, range_cont_param,
    t_tr, t_lse)

    #t_tr = 3000;
    #t_lse = 1000;

    #integ_set = (alg = Vern9(), adaptive = true, abstol = 1e-15, reltol = 1e-15);

    #params = TM6_glial_ECM_get_params();



end