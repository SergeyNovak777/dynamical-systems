username = "Alex"
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
using Pkg
Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
include(pathtorepo * "dynamical-systems\\system.jl");

using StaticArrays, JLD2, DifferentialEquations, DynamicalSystems


#sE - 1
#sI - 2
#rE -3
#rI - 4
#Y - 5
# τsE - 1, γE - 2, s0E - 3, τsI - 4, γI - 5, s0I - 6, τrE - 7, kE - 8, IE - 9,
# wEE - 10, wIE - 11, θE - 12, τrI - 13, kI - 14, II - 15,
# wEI - 16, wII - 17, θI - 18, τY - 19, βY - 20, γY - 21, ythr - 22, sEthr - 23, kY - 24

function save_output(output, pcontrolname, namefixvar, varfixvar, savevar)
    cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\brain rhythms\\diagram poincare")
    name = "poincare_map control_p $(pcontrolname) fixvar $(namefixvar) $(varfixvar) savevar $(savevar)"
    format = ".jld2"
    fullname = name * format
    jldsave(fullname; output)
end

function main()

    integ_set = (alg = RK4(), adaptive = false, dt = 0.001);

    τsE = 3.0; γE = 2.0; s0E = 0.15;
    τsI = 10.0; γI = 8.0; s0I = 0.1;

    τrE = 2.0; kE = 5.0; IE = 0.9; wEE = 3.5; wIE = 5.0; θE = 0.2;
    τrI = 6.0; kI = 5.0; II = 0.1; wEI = 5.0; wII = 3.0; θI = 0.4;

    τY = 10.0;  βY = 1.0;
    ythr = 0.5; sEthr = 0.5; kY = 0.01;
    γY = 0.0;

    p = [τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II, wEI, wII, θI, τY, βY, γY, ythr, sEthr, kY];
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0];

    ds = CoupledODEs(rate_model, u0, p, diffeq = integ_set);

    t = 100;
    ttr = 500;
    startγY = 0.0;
    endγY = 20.0;
    len = 200;
    rangeγY = range(startγY, endγY, length = len);

    name_control_param = "γY"
    index_control_param = 21;

    savevar = "rI"
    idx_save_var = 4;

    name_fix_variable = "rE";
    idx_fix = 3; fixed_value = 0.1 
    surface = (idx_fix, fixed_value)
    setting_root = (xrtol = 1e-15, atol = 1e-20)

    pmap = PoincareMap(ds, surface, rootkw = setting_root);
    output = orbitdiagram(pmap, idx_save_var, index_control_param, rangeγY; n = t, Ttr = ttr, show_progress = true);

    save_output(output, name_control_param, name_fix_variable, fixed_value, savevar);
    return 0;
end