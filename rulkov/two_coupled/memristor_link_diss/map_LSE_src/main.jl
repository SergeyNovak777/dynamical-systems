function calc_LSE(sys, type_system,
    u0, params,
    time,
    range_p1, range_p2, name_p1, name_p2,
    index_p1, index_p2)

if type_sys == "discrete"
    len_u0 = length(u0);
    prob = DiscreteProblem(sys, SVector{len_u0}(u0), time.sol, params);
    ds = DeterministicIteratedMap(sys, SVector{len_u0}(u0), params);
    len_u0 = nothing;
end

    map_LSE(prob, ds, u0, params,
    time,
    range_p1, range_p2, name_p1, name_p2,
    index_p1, index_p2)
    
end

