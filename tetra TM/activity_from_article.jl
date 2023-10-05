username = "nova"
pathtorepo = "/home/nova/work/repo_ds/dynamical-systems"
using Pkg
Pkg.activate(pathtorepo * "/env/integrate/")
#-------------------------------------------------------------------------
function TM6_glial_ECM_help()

    vars = "E - 1, x - 1, u - 2, y - 3, ECM - 4, P - 5"
    params = "τ - 1, τD - 2, τF - 3, τy - 4, α - 5, αE - 6, αecm - 7, αp - 8,\nJ - 9, U0 - 10, I0 - 11, ΔU0 - 12, β - 13, βecm  -14, βp - 15,\nγp - 16, ecm0 - 17, ecm1 - 18, kecm - 19, θecm - 20, p0 - 21, p1 - 22, θp - 23, kp - 24, ythr - 25, xthr - 26"
    return vars, params;
end

function get_params()

    τ = 0.013; τD = 0.15; τF = 1.0; τy = 1.8;   
    α = 1.5; αecm = 0.001; αp = 0.01; αE = 5.23;
    J = 3.07; U0 = 0.3; ΔU0 = 0.305; I0 = -1.741;
    β = 0.438; βp = 0.001; βecm = 0.01;
    ecm0 = 0.0; ecm1 = 1.0; kecm = 0.15; θecm = 25.6;
    p0 = 0.0; p1 = 1.0; kp = 0.05; γp = 0.1; θp = 26.0; 
    ythr = 0.5; xthr = 0.9;
    param = [τ, τD, τF, τy, α, αE, αecm, αp, J, U0, I0, ΔU0, β, βecm, βp, γp, ecm0, ecm1, kecm, θecm, p0, p1, θp, kp, ythr, xthr];
    return params;
end

function plot_timesereis(t, x, tstart, tend;
    width = 1000, height = 300, lbsize = 35, tcksize = 25,inter = true)

    f = Figure(resolution = (width, height))
    ax = Axis(f, )
    if inter
        display(GLMakie.Screen(), f)
    else
        display(f)
    end
end
#-------------------------------------------------------------------------
function main()

    init_cond = zeros(6);
    params = get_params();
    tstep = 0.001;
    integ_set = (alg = RK4(), adaptive = false, dt = tstep);


end