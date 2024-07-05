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

    if type_inheritance == "from_left_to_right" || type_inheritance == "from_down_to_up"
        preliminary_inheritance(prob, ds, Ttr, t_LSE,
            name_p1, name_p2, index_p1, index_p2, range_p1, range_p2,
            type_inheritance)
    end
end

function preliminary_inheritance(prob, ds, Ttr, t_LSE,
    name_p1, name_p2, index_p1, index_p2, range_p1, range_p2,
    type_inheritance)

    if type_inheritance == "from_left_to_right"

        name_fix_p = name_p1;
        index_fix_p = index_p1;
        value_fix_p = range_p1[1];

        name_change_p = name_p2;
        index_change_p = index_p2;
        range_change_p = range_p2;

    elseif type_inheritance == "from_down_to_up"

        name_fix_p = name_p2;
        index_fix_p = index_p2;
        value_fix_p = range_p2[2];

        name_change_p = name_p1;
        index_change_p = index_p1;
        range_change_p = range_p1;
    end

    d1_loop(prob, ds, Ttr, t_LSE,
    name_fix_p, index_fix_p, value_fix_p,
    name_change_p, index_change_p, range_change_p);
end

function d1_loop(prob, ds, Ttr, t_LSE,
    name_fix_p, index_fix_p, value_fix_p,
    name_change_p, index_change_p, range_change_p)

    params = prob.p;
    prob = remake(prob, p = params);
    prob.p[index_fix_p] = value_fix_p;
    println("$name_fix_p : $(prob.p[index_fix_p])"); flush(stdout);
    println(">>>>>>>>>>>>>>"); flush(stdout);
    for (index, value) in enumerate(range_change_p)
        prob.p[index_change_p] = value;
        println("$name_fix_p : $(prob.p[index_fix_p])"); flush(stdout);
        println("$name_change_p : $(prob.p[index_change_p])");; flush(stdout);
        println("----------------"); flush(stdout);
        # calculate trajectory
        # take last point from trajectory
        # change u0 on last point from trajectory and calculate LSE
        # inheritance
    end
end