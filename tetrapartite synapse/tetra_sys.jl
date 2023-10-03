function tetrapart_model(var, par, t)

    E, x, u, y, ecm, p  = var;
    τ, τD, τF, τy, α, αE, αecm, αp, J, U0, I0, ΔU0, β, βecm, βp, γp, ecm0, ecm1, kecm, θecm, p0, p1, θp, kp, ythr, xthr = par;

    g = log( 1 + exp( ( (J + αE * ecm) * u * x * E + I0) / α ) );
    U = U0 + ΔU0 / ( 1 + exp( -50 * ( y - ythr ) ) );
    Hy = 1 / ( 1 + exp( -20 * ( x - xthr ) ) );
    Hecm = ecm0 - (ecm0 - ecm1) / (1 + exp( -(E - θecm) / kecm ) );
    Hp =  p0 - (p0 - p1) / (1 + exp( -(E - θp) / kp) );

    dE = (-E + α * g) / τ;
    dx = (1 - x) / τD  -u * x * E;
    du = (U - u) / τF  + U * (1 - u) * E;
    dy = (-y) / τy + β * Hy;
    decm = -( αecm + γp * p ) * ecm + βecm * Hecm; 
    dp = -αp * p + βp * Hp;

    return SVector(dE, dx, du, dy, decm, dp);
end

function tetra_jac(vars, params, t)

    E_sum =  params[9] + params[6] * vars[5];

    exp_ = exp( ( E_sum * vars[3] * vars[2] * vars[1] + params[11] ) / params[5] );
    
    denomE = 1.0 + exp_;

    exp50 = exp( -50.0 * ( vars[4] - params[25] ) );
    exp20 = exp( -20.0 * ( vars[2] - params[26] ) );
    exphecm = exp( -(vars[1] - params[20]) / (params[19]) );
    exphp = exp( -(vars[1] - params[23]) / (params[24]) );

    U = params[10] + params[12] / ( 1.0 + exp( -50.0 * ( vars[4] - params[25] ) ) );
    Uder = (50.0 * params[12] * exp50) / ( 1.0 + exp50 )^2;

    Hyder = 20.0 * exp20 / (1.0 + exp20)^2;

    Hecmder = (params[17] - params[18]) * exphecm / ( params[19] * (1.0 + exphecm)^2 );
    Hpder = (params[21] - params[22]) * exphp / ( params[24] * (1.0 + exphp)^2 );

    EE = ( -1.0 + ( E_sum * vars[3] * vars[2] * exp_ ) / denomE ) / params[1];
    Ex = ( ( E_sum * vars[3] * vars[1] * exp_ ) / denomE ) / params[1];
    Eu = ( ( E_sum * vars[2] * vars[1] * exp_ ) / denomE ) / params[1];
    Eecm = ( (exp_ * params[6] * vars[3] * vars[2] * vars[1]) / denomE ) / params[1];

    xE = -vars[3] * vars[2];
    xx = - 1.0 / params[2] - vars[3] * vars[1];
    xu = -vars[2] * vars[1];


    uE = U * (1.0 - vars[3]);
    uu = -1.0 / params[3] - U * vars[1];
    uy = Uder * ( 1.0 / params[3] + (1.0 - vars[3]) * vars[1] );

    yx = params[13] * Hyder;

    ecmE = params[14] * Hecmder;
    ecmecm = -(params[7] + params[16] * vars[6]);
    ecmp = -vars[5] * params[16];

    pE = params[15] * Hpder;

    SMatrix{6,6}(EE, xE, uE, 0.0, ecmE, pE,
        Ex, xx, 0.0, yx, 0.0, 0.0,
        Eu, xu, uu, 0.0, 0.0, 0.0,
        0.0, 0.0, uy, - 1.0 / params[4], 0.0, 0.0,
        Eecm, 0.0, 0.0, 0.0, ecmecm, 0.0,
        0.0, 0.0, 0.0, 0.0, ecmp, -params[8])
end