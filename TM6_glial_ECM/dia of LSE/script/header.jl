function find_zero(Λs; tol = 1e-3)
    vec_bool = Vector{Bool}(undef, 6);
    for (index, lp) in enumerate(Λs)
        checknull=  isapprox(lp, 0.0; atol = tol);
        vec_bool[index] = checknull;
    end
    vec_of_true = findall(x->x==true, vec_bool);
    if length(vec_of_true) >=1
        return 0;
    else
        return -1;
    end;
end

function lyapunovspectrum_adapt(ds, time_lse)
    coeff_adapt_time = 1;
    flag = true;
    tands = TangentDynamicalSystem(ds; J = TM6_glial_ECM_jac);
    Λs = lyapunovspectrum(tands, time_lse);
    while flag == true
        checkzero = find_zero(Λs);
        if checkzero == -1
            coeff_adapt_time*=2;
            Λs = lyapunovspectrum(tands, time_lse * coeff_adapt_time*10);
        else
            flag = false;
        end
    end
    return Λs, coeff_adapt_time;
end