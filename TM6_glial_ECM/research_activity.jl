username = "nova";
pathtorepo = "/home/nova/work/repo_ds/dynamical-systems";
using Pkg;
Pkg.activate(pathtorepo * "/env/integrate/");

using LinearAlgebra,StaticArrays, DifferentialEquations, DynamicalSystems, IntervalRootFinding;
using CairoMakie, GLMakie;

include("/home/nova/work/repo_ds/dynamical-systems/system.jl");
include("/home/nova/work/repo_ds/dynamical-systems/visual.jl");

function TM6_glial_ECM_search_fp(ds)
    Er = interval(0.0, 30.0);
    xr, ur, yr = interval(0.0, 1.0), interval(0.0, 1.0), interval(0.0, 1.0);
    ecmr, pr = interval(0.0, 0.1), interval(0.0, 0.1);
    box = IntervalBox(Er, xr, ur, yr, ecmr, pr);
    
    fp, eigs, _ = fixedpoints(ds, box, TM6_glial_ECM_jac, tol = 1e-10,
    method = IntervalRootFinding.Krawczyk);
    return [fp, eigs]
end


function main()
time = 3000.0; tt = 0.0; tstep = 0.001; ttr = 2500.0; times = [time, tt];
integ_set = (alg = RK4(), adaptive = false, dt = tstep);

u01 = [8.8746, 0.4815, 0.8089, 0.0, 0.0, 0.0];
u02 = zeros(6);
param = TM6_glial_ECM_get_params();

ds1 = CoupledODEs(TM6_glial_ECM, u01, param, diffeq = integ_set);
tr1, trange1 = trajectory(ds1, time, Δt = tstep);

ds2 = CoupledODEs(TM6_glial_ECM, u02, param, diffeq = integ_set);
tr2, trange2 = trajectory(ds2, time, Δt = tstep);

fp, eigens = TM6_glial_ECM_search_fp(ds1);

timeseries2c = plot_timesereis_2c(trange1, trange2, tr1[:, 1], tr2[:, 1], 1, 3000000,  plot = true, width = 1600, height = 350, inter=false, lw = 1.0, color1 = :red, color2 = :blue);

data1 = [tr1[:, 2], tr1[:, 1], tr1[:, 4]];
attractor1 = plot_3d(data1, 1000000, 1500000; plot = true, prot = 60, azim = -0.55pi, elev = 0.08pi, color = :red, xl = "x", yl = "E", zl = "y", inter = false);

data2 = [tr2[:, 2], tr2[:, 1], tr2[:, 4]];
attractor2 = plot_3d(data2, 1000000, 3000000; plot = true,  prot = 60, azim = -0.55pi, elev = 0.07pi, color = :blue, xl = "x", yl = "E", zl = "y", inter = false);

idx, idy, idz = 3, 1, 4;
datamerg = [ tr1[:, idx], tr2[:, idx], tr1[:, idy], tr2[:, idy], tr1[:, idz], tr2[:, idz] ];
attractors = plot_3d_2c_fp(datamerg, [fp, idx, idy, idz], 1000000, 3000000; plot = true, prot = 60, azim = -0.55pi, elev = 0.07pi,  xl = "u", yl = "E", zl = "y", inter = false);

images = [timeseries2c, attractor1, attractor2, attractors];

return images;
end