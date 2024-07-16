include(joinpath(@__DIR__, "total_func.jl"))
#include("/home/sergey/work/repo/dynamical-systems/Map_and_diagram_LSE/total_func.jl")
#=
type_inheritance 
* from_left_to_right
* from_down_to_up
* symmetrical
* without
=#
function map_LSE_linear(sys, parameters, u0, integrator_setting,
    Ttr, t_LSE,
    name_p1, name_p2, index_p1, index_p2, range_p1, range_p2,
    type_inheritance)

    prob = ODEProblem(sys, u0, (0.0, Ttr), parameters);
    ds = CoupledODEs(sys, u0, parameters, diffeq = integrator_setting);

    dim = length(u0);
    length_p1 = length(range_p1);
    length_p2 = length(range_p2);
    
    global u0s = zeros(length_p1, length_p2, dim*2);
    global λs = zeros(length_p1, length_p2, dim);
    global namefile_LSE, namefile_u0s = get_file_name(length_p1, length_p2, name_p1, name_p2);

    dim = nothing; length_p1 = nothing; length_p2 = nothing;
    Ttr = nothing; GC.gc();

    if type_inheritance == "from_left_to_right" || type_inheritance == "from_down_to_up"
        preliminary_inheritance(prob, ds, integrator_setting,
            t_LSE,
            name_p1, name_p2, index_p1, index_p2, range_p1, range_p2,
            type_inheritance)
    end
end

function preliminary_inheritance(prob, ds, integrator_setting,
    t_LSE,
    name_p1, name_p2, index_p1, index_p2, range_p1, range_p2,
    type_inheritance)

    if type_inheritance == "from_left_to_right"

        name_fix_p = name_p1;
        index_fix_p = index_p1;
        value_fix_p = range_p1[1];

        name_change_p = name_p2;
        index_change_p = index_p2;
        range_change_p = range_p2;

        flag_save = "fix_p1";

    elseif type_inheritance == "from_down_to_up"

        name_fix_p = name_p2;
        index_fix_p = index_p2;
        value_fix_p = range_p2[2];

        name_change_p = name_p1;
        index_change_p = index_p1;
        range_change_p = range_p1;

        flag_save = "fix_p2";
    end

    d1_loop(prob, ds, integrator_setting, t_LSE,
    name_fix_p, index_fix_p, value_fix_p,
    name_change_p, index_change_p, range_change_p,
    flag_save);
end

function d1_loop(prob, ds, integrator_setting, t_LSE,
    name_fix_p, index_fix_p, value_fix_p,
    name_change_p, index_change_p, range_change_p,
    flag_save)

    if flag_save == "fix_p1"
        save_matrix_fix(change_index, u0, last_point, LSE) = save_in_matrix(1, change_index, u0, last_point, LSE)
    else
        save_matrix_fix(change_index, u0, last_point, LSE) = save_in_matrix(change_index, 1, u0, last_point, LSE)
    end
    
    last_point = prob.u0;
    params = prob.p;
    params[index_fix_p] = value_fix_p;
    
    prob = remake(prob, p = params);
    ds = reinit!(ds, p = params);

    println("$name_fix_p : $(prob.p[index_fix_p])"); flush(stdout);
    println(">>>>>>>>>>>>>>"); flush(stdout);

    for (index_cycle, value) in enumerate(range_change_p)

        prob = remake(prob, u0 = last_point);
        prob.p[index_change_p] = value;
        println("$name_fix_p : $(prob.p[index_fix_p])"); flush(stdout);
        println("$name_change_p : $(prob.p[index_change_p])");; flush(stdout);
        println("u0: $(prob.u0)"); flush(stdout);

        # calculate trajectory
        last_point = solver(prob, integrator_setting);
        
        # calculate LSE
        set_parameter!(ds, index_change_p, value);
        LSE = lyapunovspectrum(ds, t_LSE, u0 = last_point)

        # save in matrix LSE and u0s
        save_matrix_fix(index_cycle, u0, last_point, LSE)
        # inheritance
        println("----------------"); flush(stdout);
    end
end