if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
end

include("/home/sergey/work/repo/dynamical-systems/system.jl");
using JLD2, DifferentialEquations, StaticArrays;


Λs = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/R3/LSE_350x350_k_1_k_2.jld2")["λs"]
u0s = load("/home/sergey/MEGA/dynamical-systems/FHN_Korotkov/data/maps_LSE/R3/u0s_350x350_k_1_k_2.jld2")
init_point = u0s["init_points"]
last_point = u0s["last_points"]

sys = FHN2_try3;
params = FHN2_try3_params();
u0 = last_point[1, 1, :];
t_sol = 50_000;
t_span = (0.0, t_sol);

length_range = 350;
k1range = range( 0.085, 0.094, length = length_range);
k2range = range(60.0, 80.0, length = length_range);

index_parameter_1 = 7;
index_parameter_2 = 8;

name_parameter_1 = "k_1";
name_parameter_2 = "k_2";

abstl = 1e-11; reltl = 1e-11; maxiters = 10e6;
save_idxs = [1, 4];
integrator_setting = (alg = DP8(), adaptive = true, abstol = abstl, reltol = reltl, maxiters = maxiters,
save_idxs = save_idxs);
abstl = nothing; reltl = nothing; maxiters = nothing;

prob = ODEProblem(sys, SVector{length(u0)}(u0), t_span, params)


function calculate_map_EEs()

    for (index_p2_cycle, value_p2) in enumerate(k2range)
        for (index_p1_cycle, value_p1) in enumerate(k1range)

                u0 = last_point[index_p1_cycle, index_p2_cycle, :];
                prob = reinit_prob(prob, value_p1, value_p2, index_parameter_1, index_parameter_2, params, u0);
                sol = get_solve(prob, integrator_setting);
                
        end
    end

end

function reinit_prob(prob, value_p1, value_p2, index_parameter_1, index_parameter_2,
     params, u0)

    copyparams = copy(params);
    copyparams[index_parameter_1] = value_p1;
    copyparams[index_parameter_2] = value_p2;
    probcopy = remake(prob, u0 = u0, p = copyparams);
    return probcopy;
end

function get_solve(prob, integrator_setting)
    if integrator_setting.adaptive == true
        sol = solve(prob, alg = integrator_setting.alg, adaptive = true,
        abstol = integrator_setting.abstol, reltol = integrator_setting.reltol,
        maxiters = integrator_setting.maxiters,
        save_idxs = integrator_setting.save_idxs)
    else
        sol = solve(prob, alg = integrator_setting.alg, adaptive = false,
        dt = integrator_setting.dt,
        maxiters = integrator_setting.maxiters,
        save_idxs = integrator_setting.save_idxs)
    end

    return sol;    
end