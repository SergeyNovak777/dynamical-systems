function plot2d(data; type = "lines", width = 400, height = 400,
labels = ["x", "y"], title_ = "")
    f = Figure(resolution = (width, height))
    ax = Axis(f[1, 1], xlabel = labels[1], ylabel = labels[2], title = title_)
    if type == "lines"
        lines!(ax, data[1, :], data[2, :])
    elseif type == "scatter"
        scatter!(ax, data[1, :], data[2, :], markersize = 5)
    end
    display(GLMakie.Screen(), f)
end
function plot3d(data; indexs = [1, 2, 3],
width = 500, height = 500, labels = ["x", "y", "z"])
    f = Figure(resolution = (width, height))
    ax = Axis3(f[1, 1], xlabel = labels[1], ylabel = labels[2], zlabel = labels[3])
    lines!(ax, data[indexs[1], :], data[indexs[2], :], data[indexs[3], :], linewidth = 1.0)
    display(GLMakie.Screen(), f)
end
function plot3d_inter(data; indexs = [1, 2, 3],
width = 500, height = 500, scales = [1, 1, 1])
    
       f = Figure(resolution = (width, height))
       ax = ax = LScene(f[1, 1], show_axis = true)
       scale!(ax.scene, scales[1], scales[2], scales[3])
       lines!(ax, data[indexs[1], :], data[indexs[2], :], data[indexs[3], :], linewidth = 1.0)
       display(GLMakie.Screen(), f)
   end
function load_hom_curve()
    pathtofile = "C:\\Users\\" *username *  "\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab\\";
    I0 = load(pathtofile * "I0_hom_hom.jld")["data"];
    u0= load(pathtofile * "U0_hom_hom.jld")["data"];
    I0= I0[:];
    u0 = u0[:];
    return I0, u0
end
function get_norm_αs_(u, fp, A)

    count_rows, count_cols = size(u)
    index = 1
    last_index = count_rows
    matrix_αs = zeros(count_rows, count_cols)
    array_norms = zeros(count_rows)

    while index <= last_index
        norm_, αs = get_norm_αs(u[index, :], fp, A)
        matrix_αs[index, :] = αs
        array_norms[index] = norm_
        index+=1
    end
    return array_norms, matrix_αs
end

function get_norm_αs(u, fp, A)
    x = Vector(u)
    b = x - fp
    b = Vector(b)
    linprob = LinearProblem(A, b)
    linsolve = solve(linprob)
    return norm(linsolve.u, Inf), linsolve
end

function get_fp(system, jac_system, p, box,  dim = 3)
    ds = CoupledODEs(system, zeros(dim), p)
    fp, _, _  = fixedpoints(ds, box, jac_system)
    return fp
end

function get_matrix(fp, p, jac_system, norm = Inf)

    Jfp = jac_system(fp, p, 0.0)
    eigen_values_vectors = eigen(Jfp)
    eigen_vectors = eigen_values_vectors.vectors

    v1 = real(eigen_vectors[:, 1])
    v2 = real(eigen_vectors[:,2])
    v3 = imag(eigen_vectors[:, 3])

    v1 = normalize(v1, norm)
    v2 = normalize(v2, norm)
    v3 = normalize(v3, norm)

    v1 = transpose(v1)
    v2 = transpose(v2)
    v3 = transpose(v3)

    A = transpose([v1; v2; v3])

    return A
end

get_α(n, N, ϵ_box) = -ϵ_box + 2 * ϵ_box * (n-1) / (N-1)

function left_right_side(ϵ_box, count_points, fp, A, dim)
    array_points_left = zeros(count_points, dim);
    array_points_right = zeros(count_points, dim);
    n = 1;
    α1 = 0.0;
    α_vector = [α1, 0.0, 0.0];

    while n <= count_points

        α_vector[3] = get_α(n, count_points, ϵ_box);

        α_vector[2] = -ϵ_box;
        array_points_left[n, :] = fp + A * α_vector

        α_vector[2] = ϵ_box;
        array_points_right[n, :] = fp + A * α_vector
        n+=1

    end
    return array_points_left, array_points_right
end

function up_down_side(ϵ_box, count_points, fp, A, dim)
    array_points_up= zeros(count_points, dim);
    array_points_down = zeros(count_points, dim);
    n = 1;
    α1 = 0.0;
    α_vector = [α1, 0.0, 0.0]

    while n <= count_points

        α_vector[2] = get_α(n, count_points, ϵ_box);

        α_vector[3] = ϵ_box;
        array_points_up[n, :] = fp + A * α_vector

        α_vector[3] = -ϵ_box;
        array_points_down[n, :] = fp + A * α_vector
        n+=1

    end
    return array_points_up, array_points_down
end

function make_event(fp, ϵ_box, A)
    function condition(u, t, integrator)
        x = Vector(u)
        b = x - fp
        b = Vector(b)
        linprob = LinearProblem(A, b)
        linsolve = solve(linprob)
        return norm(linsolve.u, Inf) - ϵ_box
    end
    return condition
end


affect!(integrator) = terminate!(integrator)

function get_range(ϵ_shift, index_point_from_hom, p_len; param = "I0", direction_shift  = "increase")
    
    I0_hom, U0_hom = load_hom_curve();

    if param == "I0"
        pstart = I0_hom[index_point_from_hom];
    elseif param == "U0"
        pstart = I0_hom[index_point_from_hom];
    end

    if direction_shift == "increase"
        pend = pstart + ϵ_shift;
    elseif direction == "decrease"
        pend = pstart - ϵ_shift;
    end
    
    prange = range(pstart, pend, length = p_len);

    if param == "I0"
        return prange, U0_hom[index_point_from_hom];
    elseif param == "U0"
        return prange, I0_hom[index_point_from_hom];
    end
end

function get_set_p(pcontrol, pfix, pcontrolname)

    if pcontrolname == "I0"
        p = SA[α, τ, τD, τy, J, xthr, ythr, pfix, ΔU0, β, pcontrol];
    elseif pcontrolname == "U0"
        p = SA[α, τ, τD, τy, J, xthr, ythr, pcontrol, ΔU0, β, pfix];
    end

end