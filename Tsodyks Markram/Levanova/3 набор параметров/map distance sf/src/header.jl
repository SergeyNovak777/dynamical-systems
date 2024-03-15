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

using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2;

E, x, y = 0..30, 0..1, 0..1
box = E × x × y

function output(idx_I0, idx_U0, I0, U0, u0, distance)
    println("index I0: $idx_I0; I0: $I0");flush(stdout);
    println("index U0: $idx_U0; U0: $U0");flush(stdout);
    println("initial condition: $u0");flush(stdout);
    println("distance: $(distance)");flush(stdout);
end

function separate()
    println();flush(stdout);
    println("---------");flush(stdout);
end

#----------------------------------------------------

function init_ds_(params, index_control_p, control_p,index_fix_p, fix_p, u0_lc, integ_set)

    params[index_control_p] = control_p;
    params[index_fix_p] = fix_p;
    ds = CoupledODEs(TM, u0_lc, params, diffeq = integ_set);
    return ds
end

function init_ds(params, index_p1, index_p2, p1, p2, u0_lc, integ_set)
    params[index_p1] = p1;
    params[index_p2] = p2;
    ds = CoupledODEs(TM, u0_lc, params,  diffeq = integ_set)
    return ds
end

function distance_two_points_3d(p1, p2)
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    dist = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
    return dist
end

function calculate_distance(ds, t_integrate)
    
    traj, _ = trajectory(ds, t_integrate)
    fp, ei, _ = fixedpoints(ds, box, jacob_TM_);

    len_traj = length(traj)
    distance = zeros(len_traj)
    if length(fp) == 1
        for i in range(1, len_traj, step = 1)
            distance[i] = distance_two_points_3d(traj[i], fp[1])
        end
        dis = minimum(distance)
    else
        dis = -1;
    end
    return dis
end

function save_output(idx_U0, dis)
    dis_s[1, idx_U0] = dis
end

function save_output(index_I0, idx_U0, dis)
    dis_s[index_I0, idx_U0] = dis
end

function save_tofile(namefile_dis)
    jldsave(namefile_dis; dis_s)
end
