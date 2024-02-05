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

function plot_timeseries(tr, trange, title_name)
    tstart, tend = 1, 300000
    tickssize = 30
    labelsize = 40

    f= Figure(size = (1600, 400))
    axisrHz  = Axis(f[1, 1], xlabel = L"time", ylabel = L"r[Hz]", xlabelsize = labelsize, ylabelsize = labelsize,
    xticklabelsize = tickssize, yticklabelsize = tickssize, xgridvisible = false, ygridvisible = false,
    title = title_name, titlesize = labelsize)
    lines!(axisrHz, trange[tstart:tend], tr[tstart:tend, 3], linewidth = 3.0, color = :red,
    label = L"r_E")
    lines!(axisrHz, trange[tstart:tend], tr[tstart:tend, 4], linewidth = 3.0, color = :blue,
    label = L"r_I")
    axislegend(position = :rt, labelsize = 40)
    display(f)
end

function timeseries(ds, time_integrate, tstep, value_param, title_name)
    set_parameter!(ds, 21, value_param)
    tr, trange = trajectory(ds, time_integrate, Δt = tstep)
    plot_timeseries(tr, trange, title_name)
end

function timeseries_from_article()

    time_integrate = 1000.0; tstep = 0.001;
    integrate_setting = (alg = RK4(), adaptive = false, dt = 0.01)

    u0 = ones(5)*0.01;
    params = rate_model_get_params();
    value_param_array = [0.0, 25.0, 50.0, 75.0, 100.0, 125.0]
    title_name_array = [L"γY = 0.0", L"γY = 25.0", L"γY = 50.0", L"γY = 75.0",
    L"γY = 100.0", L"γY = 125.0"]

    ds = CoupledODEs(rate_model, u0, params, diffeq = integrate_setting)

    for index in range(1, step = 1, length = length(value_param_array))
       value_tmp = value_param_array[index];
       title_tmp = title_name_array[index]
       timeseries(ds, time_integrate, tstep, value_tmp, title_tmp)
    end

end