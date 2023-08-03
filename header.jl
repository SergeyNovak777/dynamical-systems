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