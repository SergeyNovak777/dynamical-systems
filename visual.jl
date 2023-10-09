function plot_timesereis(t, x, tstart, tend;
    plot = true,
    width = 1000, height = 300, ylab = "x", lbsize = 35, tcksize = 25, lw = 2.0, inter = true, grid = true, color = :black)

    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax = Axis(f[1, 1], xlabel = L"time", ylabel = ylab,
    xlabelsize = lbsize, ylabelsize = lbsize,
    xticklabelsize = tcksize, yticklabelsize = tcksize,
    xgridvisible = grid, ygridvisible = grid)
    
    lines!(ax, t[tstart:tend], x[tstart:tend], linewidth = lw, color = color);
    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end

function plot_timesereis_2c(t1, t2, x1, x2, tstart, tend;
    plot = true,
    width = 1000, height = 300, ylab = "x", lbsize = 35, tcksize = 25, lw = 2.0, inter = true, grid = true, color1 = :blue, color2 = :red)

    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax = Axis(f[1, 1], xlabel = L"time", ylabel = ylab,
    xlabelsize = lbsize, ylabelsize = lbsize,
    xticklabelsize = tcksize, yticklabelsize = tcksize,
    xgridvisible = grid, ygridvisible = grid)
    
    lines!(ax, t1[tstart:tend], x1[tstart:tend], linewidth = lw, color = color1);
    lines!(ax, t2[tstart:tend], x2[tstart:tend], linewidth = lw, color = color2);

    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end


function plot_timesereis_2f(t1, t2, x1, x2, tstart, tend;
    plot = true,
    width = 1000, height = 300, ylab = "x", lbsize = 35, tcksize = 25, lw = 2.0, inter = true, grid = true, color1 = :blue, color2 = :red)

    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax1 = Axis(f[1, 1], xlabel = L"time", ylabel = ylab,
    xlabelsize = lbsize, ylabelsize = lbsize,
    xticklabelsize = tcksize, yticklabelsize = tcksize,
    xgridvisible = grid, ygridvisible = grid)
    
    ax2 = Axis(f[2, 1], xlabel = L"time", ylabel = ylab,
    xlabelsize = lbsize, ylabelsize = lbsize,
    xticklabelsize = tcksize, yticklabelsize = tcksize,
    xgridvisible = grid, ygridvisible = grid)
    
    lines!(ax1, t1[tstart:tend], x1[tstart:tend], linewidth = lw, color = color1);
    lines!(ax2, t2[tstart:tend], x2[tstart:tend], linewidth = lw, color = color2);

    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end
    
function plot_3d(data, ts, tf;
    plot = true,
    width = 900, height = 600,
    lb_size = 30, tck_size = 25, color = :black, lw = 1.0,
    xl = "v1", yl = "v2", zl = "v3",
    azim = 1.275pi, elev = pi/8, prot = 30, disx = 40, disy = 40, disz = 40,
    grid = true, inter = false)

    x, y, z = data

    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax = Axis3(f[1, 1], azimuth = azim, elevation = elev,
                xlabel = xl, ylabel = yl, zlabel = zl,
                xlabelsize = lb_size, ylabelsize = lb_size, zlabelsize = lb_size,
                xticklabelsize = tck_size, yticklabelsize = tck_size, zticklabelsize = tck_size,
                xgridvisible = false, ygridvisible = false, zgridvisible = false,
                protrusions = prot, xlabeloffset = disx, ylabeloffset = disy, zlabeloffset = disz)

    lines!(x[ts:tf], y[ts:tf], z[ts:tf], linewidth = lw, color = color)


    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end

function plot_3d_fp(data, fixedpoint, ts, tf;
    plot = true,
    width = 900, height = 600,
    lb_size = 30, tck_size = 25, color = :black, lw = 1.0,
    xl = "v1", yl = "v2", zl = "v3",
    azim = 1.275pi, elev = pi/8, prot = 30, disx = 40, disy = 40, disz = 40,
    grid = true, inter = false)

    x, y, z = data
    fp, idx, idy, idz = fixedpoint;
    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax = Axis3(f[1, 1], azimuth = azim, elevation = elev,
                xlabel = xl, ylabel = yl, zlabel = zl,
                xlabelsize = lb_size, ylabelsize = lb_size, zlabelsize = lb_size,
                xticklabelsize = tck_size, yticklabelsize = tck_size, zticklabelsize = tck_size,
                xgridvisible = false, ygridvisible = false, zgridvisible = false,
                protrusions = prot, xlabeloffset = disx, ylabeloffset = disy, zlabeloffset = disz)

    lines!(x[ts:tf], y[ts:tf], z[ts:tf], linewidth = lw, color = color)
    for fp_ in fp
        scatter!(ax, fp_[idx], fp_[idy], fp_[idz], color = :black, markersize = 5.0);
    end

    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end

function plot_3d_fps(data, fixedpoint, ts, tf;
    plot = true,
    width = 900, height = 600,
    lb_size = 30, tck_size = 25, color = :black, lw = 1.0,
    xl = "v1", yl = "v2", zl = "v3",
    azim = 1.275pi, elev = pi/8, prot = 30, disx = 40, disy = 40, disz = 40,
    grid = true, inter = false)

    x, y, z = data
    fp, idx, idy, idz = fixedpoint;
    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax = Axis3(f[1, 1], azimuth = azim, elevation = elev,
                xlabel = xl, ylabel = yl, zlabel = zl,
                xlabelsize = lb_size, ylabelsize = lb_size, zlabelsize = lb_size,
                xticklabelsize = tck_size, yticklabelsize = tck_size, zticklabelsize = tck_size,
                xgridvisible = false, ygridvisible = false, zgridvisible = false,
                protrusions = prot, xlabeloffset = disx, ylabeloffset = disy, zlabeloffset = disz)

    lines!(x[ts:tf], y[ts:tf], z[ts:tf], linewidth = lw, color = color)
    scatter!(ax, fp[idx], fp[idy], fp[idz], color = :black, markersize = 5.0);
    

    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end

function plot_3d_2c(data, ts, tf;
    plot = true,
    width = 900, height = 600,
    lb_size = 30, tck_size = 25, color1 = :red, color2 = :blue, lw = 1.0,
    xl = "v1", yl = "v2", zl = "v3",
    azim = 1.275pi, elev = pi/8, prot = 30, disx = 40, disy = 40, disz = 40,
    grid = true, inter = false)

    x1, x2, y1, y2, z1, z2 = data
    
    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax = Axis3(f[1, 1], azimuth = azim, elevation = elev,
                xlabel = xl, ylabel = yl, zlabel = zl,
                xlabelsize = lb_size, ylabelsize = lb_size, zlabelsize = lb_size,
                xticklabelsize = tck_size, yticklabelsize = tck_size, zticklabelsize = tck_size,
                xgridvisible = false, ygridvisible = false, zgridvisible = false,
                protrusions = prot, xlabeloffset = disx, ylabeloffset = disy, zlabeloffset = disz)

    lines!(x1[ts:tf], y1[ts:tf], z1[ts:tf], linewidth = lw, color = color1)
    lines!(x2[ts:tf], y2[ts:tf], z2[ts:tf], linewidth = lw, color = color2)

    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end

function plot_3d_2c_fp(data, fixedpoint, ts, tf;
    plot = true,
    width = 900, height = 600,
    lb_size = 30, tck_size = 25, color1 = :red, color2 = :blue, lw = 1.0,
    xl = "v1", yl = "v2", zl = "v3",
    azim = 1.275pi, elev = pi/8, prot = 30, disx = 40, disy = 40, disz = 40,
    grid = true, inter = false)

    x1, x2, y1, y2, z1, z2 = data
    fp, idx, idy, idz = fixedpoint;
    CairoMakie.activate!()
    f = Figure(resolution = (width, height))
    ax = Axis3(f[1, 1], azimuth = azim, elevation = elev,
                xlabel = xl, ylabel = yl, zlabel = zl,
                xlabelsize = lb_size, ylabelsize = lb_size, zlabelsize = lb_size,
                xticklabelsize = tck_size, yticklabelsize = tck_size, zticklabelsize = tck_size,
                xgridvisible = false, ygridvisible = false, zgridvisible = false,
                protrusions = prot, xlabeloffset = disx, ylabeloffset = disy, zlabeloffset = disz)

    lines!(x1[ts:tf], y1[ts:tf], z1[ts:tf], linewidth = lw, color = color1)
    lines!(x2[ts:tf], y2[ts:tf], z2[ts:tf], linewidth = lw, color = color2)

    for fp_ in fp
        scatter!(ax, fp_[idx], fp_[idy], fp_[idz], color = :black);
    end

    if plot == true
        if inter == true
            display(GLMakie.Screen(), f)
        else
            display(f)
        end
    end
end