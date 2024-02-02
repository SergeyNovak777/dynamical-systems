if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
    using DifferentialEquations, DynamicalSystems, CairoMakie
    include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\system.jl");
else
    username = "nova"
    pathtorepo = "/home/nova/work/repo_ds/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
    using DifferentialEquations, DynamicalSystems, CairoMakie
end

integrate_setting = (alg = RK4(), adaptive = false, dt = 0.01)