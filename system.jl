#---------------------------------------------------------------------------------------------
# Tsodyks Markram
@inbounds function TM(u, p, t)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) )
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) )
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) )
    
    U_ = U(u[3], p)
    du1 = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) ) / p[2]
    du2 = (1.0 - u[2]) / p[3] - U_*u[2]*u[1]
    du3 = (-u[3])/p[4] + p[10] * σ(u[2], p)
    
    return SVector(du1, du2, du3)
end
# for bifurcationkit
@inbounds function TM_bk(u, p)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) )
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) )
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) )
    
    U_ = U(u[3], p)
    du1 = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) ) / p[2]
    du2 = (1.0 - u[2]) / p[3] - U_*u[2]*u[1]
    du3 = (-u[3])/p[4] + p[10] * σ(u[2], p)
    
    return [du1, du2, du3]
end
# for DAE
@inbounds function TM_DAE(du, u, p, t)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) )
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) )
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) )
    
    U_ = U(u[3], p)
    du[1] = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) )
    du[2] = (1.0 - u[2]) / p[3] - U_*u[2]*u[1]
    du[3] = (-u[3])/p[4] + p[10] * σ(u[2], p)
    
    nothing
end
function TM_model_get_params()
    #fixed parameters
    τ = 0.013; τD = 0.07993; τy = 3.3; J = 3.07; β = 0.300
    xthr = 0.75; ythr = 0.4
    α = 1.58; ΔU0 = 0.305;
    # control parameters
    I0 = 0.0; U0 = 0.0
    params = [α, τ, τD, τy, J, xthr, ythr, U0, ΔU0, β, I0]

    return params
end
# Jacobian
@inbounds function jacob_TM_(u, p, t)
    
    U(y, p, exp50) = p[8] + p[9] / ( 1.0 + exp50 )
    U_y(y, p, exp50) = (50.0 * p[9] * exp50) / (1.0 + exp50)^2
    g(E, x, y, p, U_) = exp((p[5]  * U_ * x * E + p[11]) / p[1])
    σ_der(x, p) = exp( (-20.0) * (x - p[6]) )
    exp50 = exp(-50.0 * (u[3] - p[7]))
    
    U_ = U(u[3], p, exp50)
    Uy = U_y(u[3], p, exp50)
    g_ = g(u[1], u[2], u[3], p, U_)
    σ_deri = σ_der(u[2], p)
    
    g_plus = 1.0 + g_
    g_mult = g_ * U_
    g_plus_mult = p[2] * (g_plus)
    u1p5 = p[5] * u[1]
    Uyu2 = Uy * u[2]
    
    E_E = (-1.0 + ((p[5] * u[2] * g_mult)) / (g_plus) ) / p[2]
    E_x = (u1p5 * g_mult) / (g_plus_mult)
    E_y = (u1p5 * Uyu2 * g_) / (g_plus_mult)
    
    x_E = -U_ * u[2]
    x_x = -1.0 / p[3] - U_ * u[1]
    x_y = -Uyu2 * u[1]
    
    y_x = 20.0 * p[10] * σ_deri / (1.0 + σ_deri)^2
    y_y = -1.0/p[4]
    
    SMatrix{3,3}(E_E, x_E, 0.0,
        E_x, x_x, y_x,
        E_y, x_y, y_y)
end

#---------------------------------------------------------------------------------------------
# Hindmarsh Rose with memristor
function HR_mem(u, p, t)
    function sigma(x)
        return 1.0 / ( 1.0 + exp( -10.0 * ( x  - ( - 0.25 ) ) ) )
    end
    memristor(z, k1_me, k2_me) = k1_me + k2_me * z^2

    a, b, c, d, s, xr, r,  I, vs, k1, k2, k1_me, k2_me  = p
    x1, y1, z1, x2, y2, z2, z = u
    
    du1 = y1 + b * x1 ^ 2 - a * x1 ^3 - z1 + I - k1 * ( x1 - vs ) * sigma(x2) + memristor(z, k1_me, k2_me)*(x2 - x1)
    du2 = c - d * x1 ^2 - y1
    du3 = r * ( s * ( x1 - xr ) - z1 )
    
    du4 = y2 + b * x2 ^ 2 - a * x2 ^3 - z2 + I - k2 * ( x2 - vs ) * sigma(x1) + memristor(z, k1_me, k2_me)*(x1 - x2)
    du5 = c - d * x2 ^2 - y2
    du6 = r * ( s * ( x2 - xr ) - z2 )

    du7 = x1 - x2
    
    return SVector(du1, du2, du3, du4, du5, du6, du7)
end
# ---------------------------------------------------------------------------------------------
# rate model
function rate_model(u, p, t)

    sE, sI, rE, rI, Y = u
    τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II, wEI, wII, θI, τY, βY, γY, ythr, sEthr, kY = p

    g(Y) = 1 + γY / (1 + exp(-Y + ythr))
    HevY(sE) = 1/(1 + exp(-(sE - sEthr)/kY))

    dsEdt = (1/τsE) *(-sE + γE*rE*(1-sE)*g(Y) + s0E)
    dsIdt = (1/τsI) *(-sI + γI*rI*(1-sI) + s0I)

    drEdt = (1/ τrE) *(-rE + 1/(1 + exp(-kE*((IE + wEE*sE-wIE*sI) -  θE))))
    drIdt = (1/τrI) *(-rI + 1/(1 + exp(-kI*((II + wEI*sE-wII*sI) - θI))))
    
    dYdt = -Y / τY + βY * HevY(sE)
    
    return SVector(dsEdt, dsIdt, drEdt, drIdt, dYdt)
end

function optim_rate_model(u, p, t)
    sE, sI, rE, rI, Y = u
    τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II, wEI, wII, θI, τY, βY, γY, ythr, sEthr, kY = p

    g(Y) = 1.0 + γY / (1.0 + exp(-Y + ythr))
    HevY(sE) = 1.0 / (1.0 + exp(-(sE - sEthr) / kY))

    dsEdt = (1.0 / τsE) * (-sE + γE * rE * (1.0 - sE) * g(Y) + s0E)
    dsIdt = (1.0 / τsI) * (-sI + γI * rI * (1.0 - sI) + s0I)

    drEdt = (1.0 / τrE) * (-rE + 1.0 / (1.0 + exp(-kE * ((IE + wEE * sE - wIE * sI) - θE))))
    drIdt = (1.0 / τrI) * (-rI + 1.0/(1.0 + exp(-kI * ((II + wEI * sE - wII * sI) - θI))))
    
    dYdt = -Y / τY + βY * HevY(sE)
    
    return SVector(dsEdt, dsIdt, drEdt, drIdt, dYdt)
end

function rate_model_get_params()
    s0E = 0.15; s0I = 0.1;
    τrE = 2.0; τrI = 6.0; τsE = 3.0; τsI = 10.0;
    IE = 0.9; II = 0.0;
    γE = 4.0; γI = 8.0;
    θE = 0.2; θI = 0.4;
    kE = 5.0; kI = 5.0;
    wEE = 3.5; wEI = 5.0; wII = 3.0; wIE = 5.0;
    γY = 0.305; βY = 0.01; τY = 0.01; kY = 0.01
    sEthr = 0.5; ythr = 0.5;

    param = [τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II, wEI, wII, θI, τY, βY, γY, ythr, sEthr, kY];
    return param
end

function rate_model_help(param)
    indexparams = "τsE = 1, γE = 2, s0E = 3, τsI = 4, γI = 5, s0I = 6, τrE = 7, kE = 8, IE = 9, wEE = 10,
    wIE = 11, θE = 12, τrI = 13, kI = 14, II = 15, wEI = 16, wII = 17, θI = 18, τY = 19, βY = 20, γY = 21, ythr = 22, sEthr = 23, kY = 24";
    nameparam = "τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II, wEI, wII, θI, τY, βY, γY, ythr, sEthr, kY";
    keyp = split(nameparam, ", ");
    dict = Dict(zip(keyp, param));

    return dict, indexparams;

end

# rate model for bifurcationkit
function rate_model(u, p)

    sE, sI, rE, rI, Y = u
    τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II, wEI, wII, θI, τY, βY, gammaY, ythr, sEthr, kY = p

    g(Y) = 1 + gammaY / (1 + exp(-Y + ythr))
    HevY(sE) = 1/(1 + exp(-(sE - sEthr)/kY))

    dsEdt = (1/τsE) *(-sE + γE*rE*(1-sE)*g(Y) + s0E)
    dsIdt = (1/τsI) *(-sI + γI*rI*(1-sI) + s0I)

    drEdt = (1/ τrE) *(-rE + 1/(1 + exp(-kE*((IE + wEE*sE-wIE*sI) -  θE))))
    drIdt = (1/τrI) *(-rI + 1/(1 + exp(-kI*((II + wEI*sE-wII*sI) - θI))))
    
    dYdt = -Y / τY + βY * HevY(sE)
    
    return [dsEdt, dsIdt, drEdt, drIdt, dYdt]
end
# jacobian rate model
function rate_jac(u, p , t)
    
    #sE, sI, rE, rI, Y = u
    τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II, wEI, wII, θI, τY, βY, γY, ythr, sEthr, kY = p

    gY = 1.0 + γY / (1.0 + exp(-u[5] + ythr))
    HevY = 1.0/(1.0 + exp(-(u[1] - sEthr)/kY))

    expE = exp(-kE*( (IE + wEE*u[1]-wIE*u[2]) - θE) )
    expI = exp( -kI*( (II + wEI*u[1]-wII*u[2]) - θI) )
    expHev = exp( -(u[1] - sEthr)/kY )
    expY = exp(-u[5] + ythr)


    sEsE = 1.0/τsE * (-1.0 - γE * u[3] *  gY)
    #sIsE = 0.0
    rEsE = 1.0/τrE * (kE * wEE * expE) / (1.0 + expE)^2
    rIsE = 1.0/τrI * (kI * wEI * expI) / (1.0 + expI)^2
    YsE = (βY * expHev) / ((1.0 + expHev)^2 * kY)

    #sEsI = 0.0
    sIsI = (1.0/τsI) * (-1.0 - γI * u[4] )
    rEsI = (1.0/τrE) * (-kE * wIE * expE) / (1.0 + expE)^2 
    rIsI = (1.0/τrI) * (-kI * wII * expI) / (1.0 + expI)^2
    #YsI = 0.0

    sErE = (1.0/τsE) * ( γE*(1-u[1])*gY )
    #sIrE = 0.0
    rErE = (-1.0/τrE)
    #rIrE = 0.0
    #YrE = 0.0

    #sErI = 0.0
    sIrI = (1.0/τsI) * ( γI*(1.0-u[2]) )
    #rErI = 0.0
    rIrI = (-1.0/τrI)
    #YrI = 0.0

    sEY = (1.0/τsE) * ( γE*u[3]*(1.0-u[1]) * ( γY * expY ) / (1.0 + expY )^2 )
    #sIY = 0.0
    #rEY = 0.0
    #rIY = 0.0
    YY = (-1.0/τY)

    return SMatrix{5,5}(sEsE, 0.0, rEsE, rIsE, YsE,
    0.0, sIsI, rEsI, rIsI, 0.0,
    sErE, 0.0, rErE, 0.0, 0.0,
    0.0, sIrI, 0.0, rIrI, 0.0,
    sEY, 0.0, 0.0, 0.0, YY)

end
#---------------------------------------------------------------------------------

function arneodo_model(u, p, t)
    x = u[2];
    y = u[3];
    z = -u[2] - p[1] * u[3] + p[2] * u[1] * ( 1 - u[1] );
    return SVector(x, y ,z);
end

function jac_arneodo(u, p, t)
    return SMatrix{3, 3}(
    0.0, 0.0, p[2]*(1-2*u[1]),
    1.0, 0.0, -1.0,
    0.0, 1.0, -p[1];)
end

# ---------------------------------------------------------------------------------
# model tetrapartite synapse

@inbounds function TM6_glial_ECM(var, par, t)

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
    E_mult_x_mult_u = vars[3] * vars[2] * vars[1]

    exp_ = exp( ( E_sum * E_mult_x_mult_u + params[11] ) / params[5] );
    
    denomE = 1.0 + exp_;
    y_minus_ythr = vars[4] - params[25];
    E_sum_mult_exp_ = E_sum * exp_;
    E_sum_mult_exp_mult_E = E_sum_mult_exp_ * vars[1];
    

    exp50 = exp( -50.0 * ( y_minus_ythr ) );
    exp20 = exp( -20.0 * ( vars[2] - params[26] ) );
    exphecm = exp( -(vars[1] - params[20]) / (params[19]) );
    exphp = exp( -(vars[1] - params[23]) / (params[24]) );

    U = params[10] + params[12] / ( 1.0 + exp( -50.0 * ( y_minus_ythr ) ) );
    Uder = (50.0 * params[12] * exp50) / ( 1.0 + exp50 )^2;

    Hyder = 20.0 * exp20 / (1.0 + exp20)^2;

    Hecmder = (params[17] - params[18]) * exphecm / ( params[19] * (1.0 + exphecm)^2 );
    Hpder = (params[21] - params[22]) * exphp / ( params[24] * (1.0 + exphp)^2 );


    EE = ( -1.0 + ( E_sum_mult_exp_ * vars[3] * vars[2] ) / denomE ) / params[1];
    Ex = ( ( E_sum_mult_exp_mult_E * vars[3]  ) / denomE ) / params[1];
    Eu = ( ( E_sum_mult_exp_mult_E * vars[2]  ) / denomE ) / params[1];
    Eecm = ( (exp_ * params[6] * E_mult_x_mult_u) / denomE ) / params[1];


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

#=@inbounds function TM6_glial_ECM_jac(vars, params, t)

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
end=#

function TM6_glial_ECM_get_params()
    τ = 0.013; τD = 0.15; τF = 1.0; τy = 1.8;   
    α = 1.5; αecm = 0.001; αp = 0.01;
    J = 3.07; U0 = 0.3; ΔU0 = 0.305; 
    β = 0.438; βp = 0.01; βecm = 0.01;
    ecm0 = 0.0; ecm1 = 1.0; kecm = 0.15; θecm = 25.6;
    p0 = 0.0; p1 = 1.0; kp = 0.05; γp = 0.1; θp = 26.0; 
    ythr = 0.5; xthr = 0.9;
    αE = 5.0;
    I0 = -1.741;

    param = [τ, τD, τF, τy, α, αE, αecm, αp, J, U0, I0, ΔU0, β, βecm, βp, γp, ecm0, ecm1, kecm, θecm, p0, p1, θp, kp, ythr, xthr];
    return param;
end
function TM6_glial_ECM_help(params)
    
    indexparams = "
    τ - 1, τD - 2, τF - 3, τy - 4, α - 5, αE - 6, αecm - 7, αp - 8,
    J - 9, U0 - 10, I0 - 11, ΔU0 - 12, β - 13, βecm  -14, βp - 15,
    γp - 16, ecm0 - 17, ecm1 - 18, kecm - 19, θecm - 20, p0 - 21, p1 - 22, θp - 23, kp - 24, ythr - 25, xthr - 26";

    nameparams = "τ, τD, τF, τy, α, αE, αecm, αp, J, U0, I0, ΔU0, β, βecm, βp, γp, ecm0, ecm1, kecm, θecm, p0, p1, θp, kp, ythr, xthr";
    keyp = split(nameparams, ", ");
    dict = Dict(zip(keyp, params));
    
    return dict, indexparams;
end

# ---------------------------------------------------------------------------------
# model two coupled FHN memristor
function FHN2_try3(u, p ,t)
    x1, y1, x2, y2, z= u
    ϵ, a, g, k, σ, α, k1, k2 = p

    I(ϕ_i) = g * (1.0/(1.0 + exp(k*(cos(σ/2) - cos(ϕ_i - α - σ/2)))))
    ρz = k1 + k2 * z ^ 2

    ϕ2 = atan(y2, x2)
    ϕ1 = atan(y1, x1)

    dx1dt = (x1 - x1 ^ 3 / 3 - y1 + I(ϕ2) + ρz * (x2 - x1) ) / ϵ
    dy1dt = x1 - a
    dx2dt = (x2 - x2 ^ 3 / 3 - y2 + I(ϕ1) + ρz * (x1 - x2) ) / ϵ
    dy2dt = x2 - a
    dzdt = x1 - x2
    return SVector(dx1dt, dy1dt, dx2dt, dy2dt, dzdt)
end

function FHN2_try3_params()
    ϵ = 0.01; a = -1.01;
    g = 0.1; k = 50.0; σ = 50.0 * pi / 180; α = 160.0 * pi / 180;
    k1 = 0.0; k2 = 0.0
    return [ ϵ, a, g, k, σ, α, k1, k2]
end