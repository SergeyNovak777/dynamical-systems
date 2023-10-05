# linux
username = "nova"
pathtorepo = "/home/nova/work/repo_ds/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")

using StaticArrays, DifferentialEquations, DynamicalSystems, ForwardDiff, BenchmarkTools, IntervalRootFinding 
using CairoMakie, GLMakie


function not_optim(var, par, t)

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

@inbounds function optim(var, par, t)

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

function main()

    time = 10000.0; tt = 0.0; tstep = 0.001; ttr = 2500.0; times = [time, tt];
    integ_set = (alg = RK4(), adaptive = false, dt = tstep);

    τ = 0.013; τD = 0.15; τF = 1.0; τy = 1.8;   
    α = 1.5; αecm = 0.001; αp = 0.01;
    J = 3.07; U0 = 0.3; ΔU0 = 0.305; 
    β = 0.438; βp = 0.001; βecm = 0.01;
    ecm0 = 0.0; ecm1 = 1.0; kecm = 0.15; θecm = 25.6;
    p0 = 0.0; p1 = 1.0; kp = 0.05; γp = 0.1; θp = 26.0; 
    ythr = 0.5; xthr = 0.9;
    αE = 5.23;
    I0 = -1.741;

    u0 = zeros(6);
    param = [τ, τD, τF, τy, α, αE, αecm, αp, J, U0, I0, ΔU0, β, βecm, βp, γp, ecm0, ecm1, kecm, θecm, p0, p1, θp, kp, ythr, xthr];
end