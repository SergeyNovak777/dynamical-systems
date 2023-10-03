function tetra_jac(vars, params, t)

    E, x, u, y, ecm, p  = vars;
    # E - 1, x - 2, u - 3, y - 4, ecm - 5, p -6
    τ, τD, τF, τy, α, αE, αecm, αp,
    J, U0, I0, ΔU0, β, βecm, βp, γp,
    ecm0, ecm1, kecm, θecm, p0, p1, θp, kp, ythr, xthr = params;

    E_sum =  J + αE * ecm;

    exp_ = exp( ( E_sum * u * x * E + I0 ) / α );
    
    denomE = 1.0 + exp_;

    exp50 = exp( -50.0 * ( y - ythr ) );
    exp20 = exp( -20.0 * ( x - xthr ) );
    exphecm = exp( -(E - θecm) / (kecm) );
    exphp = exp( -(E - θp) / (kp) );

    U = U0 + ΔU0 / ( 1.0 + exp( -50.0 * ( y - ythr ) ) );
    Uder = (50.0 * ΔU0 * exp50) / ( 1.0 + exp50 )^2;

    Hyder = 20.0 * exp20 / (1.0 + exp20)^2;

    Hecmder = (ecm0 - ecm1) * exphecm / ( kecm * (1.0 + exphecm)^2 );
    Hpder = (p0 - p1) * exphp / ( kp * (1.0 + exphp)^2 );

    EE = ( -1.0 + ( E_sum * u * x * exp_ ) / denomE ) / τ;
    Ex = ( ( E_sum * u * E * exp_ ) / denomE ) / τ;
    Eu = ( ( E_sum * x * E * exp_ ) / denomE ) / τ;
    Eecm = ( (exp_ * αE * u * x * E) / denomE ) / τ;

    xE = -u * x;
    xx = - 1.0 / τD - u * E;
    xu = -x * E;


    uE = U * (1.0 - u);
    uu = -1.0 / τF - U * E;
    uy = Uder * ( 1.0 / τF + (1.0 - u) * E );


    yy = - 1.0 / τy;
    yx = β * Hyder;

    ecmE = βecm * Hecmder;
    ecmecm = -(αecm + γp * p);
    ecmp = -ecm * γp;

    pE = βp * Hpder;
    pp = -αp;

    SMatrix{6,6}(EE, xE, uE, 0.0, ecmE, pE,
        Ex, xx, 0.0, yx, 0.0, 0.0,
        Eu, xu, uu, 0.0, 0.0, 0.0,
        0.0, 0.0, uy, yy, 0.0, 0.0,
        Eecm, 0.0, 0.0, 0.0, ecmecm, 0.0,
        0.0, 0.0, 0.0, 0.0, ecmp, pp)
end