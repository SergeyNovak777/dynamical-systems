if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
    using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2
        include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\system.jl");
else
    username = "nova"
    pathtorepo = "/home/nova/work/repo_ds/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
    using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2
    include("/home/nova/work/repo_ds/dynamical-systems/system.jl")
end

function save_logs(file_logs,
    index, name_p, p_range, p_value,
    u0_st, u0_ed)

    open(file_logs, "a") do file
        print(file, "index: $(index); $(name_p) $(p_range[index]) \n");
        print(file, "$(name_p) from ds $(name_p): $(p_value); \n");
        print(file, "u0 start: $(u0_st); \n");
        print(file, "u0 end: $(u0_ed) \n");
        print(file, "\n");
    end
end

function save_logs(file_logs,
    index_p1, index_p2, name_p1, name_p2,
    range_p1, range_p2,
    p1_value, p2_value,
    u0_st, u0_ed)

    open(file_logs, "a") do file
        print(file, "index p1 : $(index_p1); $(name_p1): $(range_p1[index_p1]) \n");
        print(file, "index p2 : $(index_p2); $(name_p2): $(range_p2[index_p2]) \n");
        print(file, "$(name_p1) from ds $(name_p1): $(p1_value); \n");
        print(file, "$(name_p2) from ds $(name_p2): $(p2_value); \n");
        print(file, "u0 start: $(u0_st); \n");
        print(file, "u0 end: $(u0_ed) \n");
        print(file, "\n");
    end
end
function get_ds()

    tol = 1e-15;
    integrator_setting = (alg = Vern9(), adaptive = true, abstol = tol, reltol = tol);
    params = TM6_glial_ECM_get_params();
    u0 = [0.9445509341100914, 0.74116702856987, 0.7361196042973006, 0.0646914552140727, 0.15145764079879162, 0.0009327645775731449];
    
    ds= CoupledODEs(TM6_glial_ECM, u0, params, diffeq = integrator_setting);
    return ds;
end

function get_control_param()
    name_p1 = "I0";
    name_p2 = "αE";
    index_p1 = 11;
    index_p2 = 6;
    return [name_p1, name_p2, index_p1, index_p2];
end

function map2d_u0()
    
    t_transient = 3000;
    global u0_loc_2d = zeros(6);

    name_p1, name_p2, index_p1, index_p2 = get_control_param();

    len = 250;
    range_p1 = range(-1.741, -1.6, length = len);
    range_p2 = range(0.067, 5.0, length = len);

    ds = get_ds();
    systems = [deepcopy(ds) for _ in 1:Threads.nthreads()-1];
    pushfirst!(systems, ds);
    
    total_file_names = "map u0s $(name_p1) $(name_p2) $(length(range_p1))*$(length(range_p2))";
    file_logs_p2 = total_file_names * " change p2 logs.txt";
    file_logs_2d = total_file_names * " change p1 p2 logs.txt";
    u0s = zeros(length(range_p1), length(range_p2), dimension(ds));

    

    for index in eachindex(range_p2)

        if index == 1
            global u0_local_prep = initial_state(ds);
        end

        u0_st = u0_local_prep;

        set_parameter!(ds, index_p2, range_p2[index]);
        tr, _ = trajectory(ds, t_transient, u0_local_prep);

        u0_local_prep = tr[end];
        u0s[1, index, :] = u0_local_prep;

        save_logs(file_logs_p2, index,
        name_p2, range_p2, current_parameters(ds)[index_p2],
        u0_st, u0_local_prep);
    end

    Threads.@threads for index_ext in eachindex(range_p2)
        for index_int in eachindex(range_p1)

            if index_int == 1
                continue;
            end

            system = systems[Threads.threadid()]

            u0_loc_2d = u0s[index_int - 1, index_ext, :];
            u0_loc_2d_st = u0_loc_2d;

            set_parameter!(system, index_p2, range_p2[index_ext]);
            set_parameter!(system, index_p1, range_p1[index_int]);
            tr, _= trajectory(system, t_transient, u0_loc_2d);

            u0_loc_2d = tr[end];

            u0s[index_int, index_ext, :] = u0_loc_2d;
            
            save_logs(file_logs_2d,
            index_int, index_ext, name_p1, name_p2,
            range_p1, range_p2,
            current_parameters(system)[index_p1], current_parameters(system)[index_p2],
            u0_loc_2d_st, u0_loc_2d);

        end
    end

end