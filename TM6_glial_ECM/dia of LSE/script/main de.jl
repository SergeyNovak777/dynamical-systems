if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
    using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2
    include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\system.jl");
    include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\TM6_glial_ECM\\dia of LSE\\script\\test_inheritance.jl");
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
        print(file, "index: $(index) $(name_p) $(p_range[index]) \n");
        print(file, "from ds $(name_p): $(p_value) \n");
        print(file, "u0 start: $(u0_st) \n");
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
        print(file, "index p1 : $(index_p1) $(name_p1): $(range_p1[index_p1]) \n");
        print(file, "index p2 : $(index_p2) $(name_p2): $(range_p2[index_p2]) \n");
        print(file, "$(name_p1) from ds $(name_p1): $(p1_value) \n");
        print(file, "$(name_p2) from ds $(name_p2): $(p2_value) \n");
        print(file, "u0 start: $(u0_st) \n");
        print(file, "u0 end: $(u0_ed) \n");
        print(file, "\n");
    end
end

function solve_(sys)

    alg = Vern9();
    adapt = true; tol = 1e-12; abstol = tol; reltol = tol;

    sol = solve(sys, alg = alg,
    adaptive = adapt, abstol = abstol, reltol = reltol);
    return sol;
end

function get_control_param()
    name_p1 = "I0";
    name_p2 = "αE";
    index_p1 = 11;
    index_p2 = 6;
    return [name_p1, name_p2, index_p1, index_p2];
end

function map2d_u0()

    u0start = SA[0.9445509341100914, 0.74116702856987, 0.7361196042973006, 0.0646914552140727, 0.15145764079879162, 0.0009327645775731449];

    D = length(u0start)
    exitcode = 0;
    u0_local_prep = zeros(D);
    u0_loc_2d = zeros(D);

    tstart = 0.0; tend = 1500.0; tspan = (tstart, tend);
    name_p1, name_p2, index_p1, index_p2 = get_control_param();

    len = 50;
    range_p1 = range(-1.741, -1.6, length = len);
    range_p2 = range(0.067, 5.0, length = len);
    
    filename_u0ttrs  = "map u0ttrs threads $(name_p1) $(name_p2) $(length(range_p1))x$(length(range_p2))";
    filename_u0s  = "map u0s threads $(name_p1) $(name_p2) $(length(range_p1))x$(length(range_p2))";
    
    #file_logs_p2 = total_file_names * " change p2 logs.txt";
    #file_logs_2d = total_file_names * " change p1 p2 logs.txt";

    u0ttrs = zeros(length(range_p1), length(range_p2), length(u0start));
    u0s = zeros(length(range_p1), length(range_p2), length(u0start));
    
    # preparation inheritance
    for index in eachindex(range_p2)

        if index == 1
            u0_local_prep = u0start;
        end

        u0_st = u0_local_prep;

        params = TM6_glial_ECM_get_params();
        params[index_p2] = range_p2[index];
        
        system = ODEProblem(TM6_glial_ECM, u0_st, tspan, params);
        sol = solve_(system);

        if sol.retcode == ReturnCode.MaxIters
            exitcode = 1;
        end

        u0_local_prep = sol[end];

        u0ttrs[1, index, :] = u0_st;
        u0s[1, index, :] = u0_local_prep;
        
        #save_logs(file_logs_p2, index,
        #name_p2, range_p2, system.p[index_p2],
        #u0_st, u0_local_prep);

        if exitcode == 1
            exit();
        end;
    end;
    jldsave(filename_u0s*".jld2"; u0s);
    jldsave(filename_u0ttrs*".jld2"; u0ttrs);

    #inheritance in map
    Threads.@threads :static for index_ext in eachindex(range_p2)
        for index_int in eachindex(range_p1)

            if index_int == 1
                u0_loc_2d = u0s[index_int, index_ext, :];
                u0_loc_2d = SVector{D}(u0_loc_2d);
                continue;
            end

            u0_loc_2d_st = u0_loc_2d;

            params = TM6_glial_ECM_get_params();

            params[index_p1] = range_p1[index_int];
            params[index_p2] = range_p2[index_ext];
            
            system = ODEProblem(TM6_glial_ECM, u0_loc_2d_st, tspan, params);
            sol = solve_(system);

            if sol.retcode == ReturnCode.MaxIters
                exitcode = 1;
            end

            u0_loc_2d = sol[end];

            u0ttrs[index_int, index_ext, :] = u0_loc_2d_st;
            u0s[index_int, index_ext, :] = u0_loc_2d;
            
            #save_logs(file_logs_2d, index_int, index_ext, name_p1, name_p2,  range_p1, range_p2, system.p[index_p1], system.p[index_p2], u0_loc_2d_st, u0_loc_2d);

            if exitcode == 1
                exit();
            end;
        end
    end

    jldsave(filename_u0s*".jld2"; u0s);
    jldsave(filename_u0ttrs*".jld2"; u0ttrs);

    checks(u0s, u0ttrs)
    return u0s, u0ttrs;
end 