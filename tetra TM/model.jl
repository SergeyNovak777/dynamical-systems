@inbounds function TM6_glial_ECM(var, params, t)

    g = log( 1.0 + exp( ( (par[9] + par[6] * var[5]) * var[3] * var[2] * var[1] + par[11]) / par[5] ) );
    U = par[10] + par[12] / ( 1.0 + exp( -50.0 * ( var[4] - par[25] ) ) );
    Hy = 1.0 / ( 1.0 + exp( -20.0 * ( var[2] - par[26] ) ) );
    Hecm = par[17] - (par[17] - par[18]) / (1.0 + exp( -(var[1] - par[20]) / par[19] ) );
    Hp =  par[21] - (par[21] - par[22]) / (1.0 + exp( -(var[1] - par[23]) / par[24]) );

    dE = (-var[1] + par[5] * g) / par[1];
    dx = (1.0 - var[2]) / par[2]  -var[3] * var[2] * var[1];
    du = (U - var[3]) / par[3]  + U * (1.0 - var[3]) * var[1];
    dy = (-var[4]) / par[4] + par[13] * Hy;
    decm = -( par[7] + par[16] * var[6] ) * var[5] + par[14] * Hecm; 
    dp = -par[8] * var[6] + par[15] * Hp;

    return SVector(dE, dx, du, dy, decm, dp);
end

@inbounds function TM6_glial_ECM_jac(vars, params, t)

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

    ecmE =  -params[14] * Hecmder;
    ecmecm = -(params[7] + params[16] * vars[6]);
    ecmp = -vars[5] * params[16];

    pE = -params[15] * Hpder;

    SMatrix{6,6}(EE, xE, uE, 0.0, ecmE, pE,
        Ex, xx, 0.0, yx, 0.0, 0.0,
        Eu, xu, uu, 0.0, 0.0, 0.0,
        0.0, 0.0, uy, - 1.0 / params[4], 0.0, 0.0,
        Eecm, 0.0, 0.0, 0.0, ecmecm, 0.0,
        0.0, 0.0, 0.0, 0.0, ecmp, -params[8])
end