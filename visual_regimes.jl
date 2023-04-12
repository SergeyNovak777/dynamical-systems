function research_regime(model, u0, p, times, integ_setting;
        index_3d = [3, 2, 1], index_ts = 1,
        t_int_ts = nothing, t_int_3d = nothing,
        plot = true, plot_phase_space = true, plot_t_series = true,
        return_point_last = false)
    
    ds = CoupledODEs(model, u0, p, diffeq = integ_setting)
    tr, trange = trajectory(ds, times[1]; Δt = integ_setting.dt, Ttr = times[2])
    
    t_int_ts = check_t_int(t_int_ts, times)
    t_int_3d = check_t_int(t_int_3d, times)

    phase_space = interactive_3d(tr, index_3d, t_int_3d)
    tseries = timeseries([trange, tr], index_ts, t_int_ts)
    
    if plot == true
        GLMakie.activate!()
        if plot_phase_space == true && plot_t_series == true
            display(GLMakie.Screen(), phase_space)
            display(GLMakie.Screen(), tseries)
        elseif plot_phase_space == true
            display(GLMakie.Screen(), phase_space)
        else
            display(GLMakie.Screen(), tseries)
        end
    end

    if return_point_last ==  false
        return [ds, tr, trange, phase_space, tseries]
    else
        return [ds, tr[end], phase_space, tseries]
    end
end

function interactive_3d(data, index, interval)
    f = Figure(resolution = (800, 500))
    axis3 = Axis3(f[1, 1],
                xlabelsize = 35, ylabelsize = 35, zlabelsize = 35,
                xticklabelsize = 20, yticklabelsize = 20, zticklabelsize = 20,
                xgridvisible = false, ygridvisible = false, zgridvisible = false)

    lines!(axis3,
        data[interval[1]:interval[2], index[1]], data[interval[1]:interval[2], index[2]],
        data[interval[1]:interval[2],index[3]], linewidth = 1.0, color = :deeppink)
    return f
end

function timeseries(data, index, interval)
    f = Figure(resolution = (900, 250))
    axis = Axis(f[1, 1],
            xlabelsize = 35, ylabelsize = 35,
            xticklabelsize = 28, yticklabelsize = 28,
            xgridvisible = false, ygridvisible = false)
    
    lines!(axis, data[1][interval[1]:interval[2]], data[2][interval[1]:interval[2], index],
            linewidth = 1.0, color = :deeppink)
    return f
end

function research_regime_2_traj(model, u0_1, u0_2, p, times, integ_setting;
    index_3d = [3, 2, 1], t_int_3d = nothing,
    plot = true, onlylastpoint = false, dis_hc = false)

    ds_1 = CoupledODEs(model, u0_1, p, diffeq = integ_setting) 
    ds_2 = CoupledODEs(model, u0_2, p, diffeq = integ_setting) 

    tr_1, trange = trajectory(ds_1, times[1]; Δt = integ_setting.dt, Ttr = times[2])
    tr_2, _ = trajectory(ds_2, times[1]; Δt = integ_setting.dt, Ttr = times[2])

    fp, eigs, _ = fixedpoints(ds_1, box, jacob_TM_)

    if dis_hc == true
        distance_1 = zeros(length(tr_1))
        distance_2 = zeros(length(tr_2))

        for i in range(1, length(distance_1), step = 1)
            distance_1[i] = distance_btw_plot(tr_1[i], fp[1])
        end

        for i in range(1, length(distance_2), step = 1)
            distance_2[i] = distance_btw_plot(tr_2[i], fp[1])
        end

        dis_1 = minimum(distance_1)
        dis_2 = minimum(distance_2)
    end

    t_int_3d = check_t_int(t_int_3d, times)

    inter_3d = interactive_3d_2_traj(tr_1, tr_2, fp[1], index_3d, t_int_3d)

    if plot == true
        GLMakie.activate!()
        display(GLMakie.Screen(), inter_3d)
    end

    if onlylastpoint == true && dis_hc == true
        return [ds_1, ds_2, tr_1[end], tr_2[end], dis_1, dis_2, inter_3d, fp, eigs]

    else
        return [ds_1, ds_2, tr_1, tr_2, dis_1, dis_2, inter_3d, fp, eigs]
    end

end

function interactive_3d_2_traj(data_1, data_2, fp, index, interval)

    f = Figure(resolution = (800, 500))
    axis3 = Axis3(f[1, 1],
            xlabelsize = 35, ylabelsize = 35, zlabelsize = 35,
            xticklabelsize = 20, yticklabelsize = 20, zticklabelsize = 20,
            xgridvisible = false, ygridvisible = false, zgridvisible = false)

    lines!(axis3,
        data_1[interval[1]:interval[2], index[1]], data_1[interval[1]:interval[2], index[2]],
        data_1[interval[1]:interval[2],index[3]], linewidth = 1.0, color = :deeppink)

    lines!(axis3,
        data_2[interval[1]:interval[2], index[1]], data_2[interval[1]:interval[2], index[2]],
        data_2[interval[1]:interval[2],index[3]], linewidth = 1.0, color = :black)

    scatter!(axis3, fp[index[1]], fp[index[2]], fp[index[3]], linewidth = 5.0, color = :orange)
    return f
end

function check_t_int(t_int, times)
    if t_int == nothing
        t_int = [1, times[1]]
    end
    return t_int
end

function distance_btw_plot(p1, p2)
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    dist = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )
    return dist
end