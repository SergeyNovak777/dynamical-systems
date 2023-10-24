if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
    using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2
    #include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\system.jl");
    #include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\TM6_glial_ECM\\dia of LSE\\script\\test_inheritance.jl");
else
    username = "nova";
    pathtorepo = "/home/nova/work/repo_ds/dynamical-systems";
    using Pkg;
    Pkg.activate(pathtorepo * "/env/integrate/");
    using StaticArrays, DifferentialEquations, DynamicalSystems, JLD2;
    include("/home/nova/work/repo_ds/dynamical-systems/system.jl");
    include("/home/nova/work/repo_ds/dynamical-systems/TM6_glial_ECM/dia of LSE/inheritance/includes/linear_map_u0s.jl");
    include("/home/nova/work/repo_ds/dynamical-systems/TM6_glial_ECM/dia of LSE/inheritance/includes/test.jl");
end

function get_control_param()
    name_p1 = "I0";
    name_p2 = "αE";
    index_p1 = 11;
    index_p2 = 6;
    return [name_p1, name_p2, index_p1, index_p2];
    # name_p1, name_p2, index_p1, index_p2
end


#--------------------------------------------------------------------
sys = TM6_glial_ECM;
u0 = SA[0.9445509341100914, 0.74116702856987, 0.7361196042973006, 0.0646914552140727, 0.15145764079879162, 0.0009327645775731449];
params = TM6_glial_ECM_get_params();
tspan = (0.0, 1500.0);
    
alg = Vern9(); adaptive = true;
tol = 1e-12; abstol = tol; reltol = tol;
integ_set = (alg = alg, adaptive = adaptive, abstol = abstol, reltol = reltol);

control_params = get_control_param();
    
len = 50;
I0_range = range(-1.741, -1.6, length = len);
αE_range = range(0.067, 5.0, length = len);

#matrix_u0_start, matrix_u0_end = 
@time linear_2pmap(sys, params, u0, tspan, integ_set, control_params, I0_range, αE_range);

check_inheritance_vertical(matrix_u0_start, matrix_u0_end, len)