username = "admin";
env = "integrate";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");
include(pathtorepo * "dynamical-systems\\system.jl");

using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2;

#function output(idx_I0, idx_U0, I0, U0, u0)
    #println("index I0: $idx_I0; I0: $I0");flush(stdout);
    #println("index U0: $idx_U0; U0: $U0");flush(stdout);
    #println("initial condition: $u0");flush(stdout);
#end

#function output(idx_U0,U0, u0)
    #println("index U0: $idx_U0; U0: $U0");flush(stdout);
    #println("initial condition: $u0");flush(stdout);
#end

function output(p1_name, p2_name, index_p1, index_p2, p1, p2, u0)
    println("index $(p1_name): $(index_p1); $(p1_name) value: $(p1)");flush(stdout);
    println("index $(p2_name): $(index_p2); $(p2_name) value: $(p2)");flush(stdout);
    println("initial condition: $u0");flush(stdout);
end

function output(p1_name, index_p1, p1, u0)
    println("index $(p1_name): $(index_p1); $(p1_name) value: $(p1)");flush(stdout);
    println("initial condition: $u0");flush(stdout);
end

function output_end_iter(ΛΛ, u0_lc)
    println("LSE: $(ΛΛ)");flush(stdout);
    println("last point: $(u0_lc)");flush(stdout);
end

function separate()
    println();flush(stdout);
    println("---------");flush(stdout);
end

#----------------------------------------------------

function init_ds_(sys, params, index_control_p, control_p,index_fix_p, fix_p, u0_lc, integ_set)

    params[index_control_p] = control_p;
    params[index_fix_p] = fix_p;
    ds = CoupledODEs(sys, u0_lc, params, diffeq = integ_set);
    return ds
end

function init_ds(sys, params, index_p1, index_p2, p1, p2, u0_lc, integ_set)
    params[index_p1] = p1;
    params[index_p2] = p2;
    ds = CoupledODEs(sys, u0_lc, params,  diffeq = integ_set)
    return ds
end

function goto_attractor(ds_, time_attract, integ_set)
    tr,_ = trajectory(ds_, time_attract; Δt = integ_set.dt)
    u0_lc = tr[end]
    return u0_lc
end

function spectrum(ds_, t)
    ΛΛ = lyapunovspectrum(ds_, t)
    return ΛΛ
end

function save_output(index_p2, ΛΛ, u0_lc)
    Λs[1, index_p2, :] = ΛΛ
    u0s[1, index_p2, :] = u0_lc
end

function save_output(index_p1, index_p2, ΛΛ, u0_lc)
    Λs[index_p1, index_p2, :] = ΛΛ
    u0s[index_p1, index_p2, :] = u0_lc
end

function save_tofile(namefile_LSE, namefile_u0s)
    jldsave(namefile_LSE; Λs)
    jldsave(namefile_u0s; u0s)
end
