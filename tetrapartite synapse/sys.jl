function tetrapart_model(var, par, t)

    g = log( 1 + exp( ( (par[9] + par[6] * var[5]) * var[3] * var[2] * var[1] + par[11]) / par[5] ) );
    U = par[10] + par[12] / ( 1.0 + exp( -50.0 * ( var[4] - par[25] ) ) );
    Hy = 1 / ( 1 + exp( -20.0 * ( var[2] - par[26] ) ) );
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