function load_params()
    τ = 0.013; τD = 0.07993; τy = 3.3; J = 3.07; β = 0.300;
    xthr = 0.75; ythr = 0.4; α = 1.58; ΔU0 = 0.305;
    p = [α, τ, τD, τy, J, xthr, ythr, 0.0, ΔU0, β, 0.0];
    M = Diagonal([τ, 1.0, 1.0])
    return [p, M]
end

function load_hom_curve()
    cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab");
    I0_hom = load("I0_hom_hom.jld")["data"];
    u0_hom = load("U0_hom_hom.jld")["data"];
    I0_hom = I0_hom[:];
    U0_hom = u0_hom[:];
    return I0_hom, U0_hom;
end

function get_fp(system, jac_system, p, box,  dim = 3)
    ds = CoupledODEs(system, zeros(dim), p)
    fp, _, _  = fixedpoints(ds, box, jac_system)
    return fp
end

function check_fp(fp)
    if length(fp) == 1
        return fp[1];
    else
        throw(DomainError(fp, "more one fixed point"));
    end
end

function get_matrix(fp, p, jac_system, norm = 2)

    Jfp = jac_system(fp, p, 0.0);
    eigen_values_vectors = eigen(Jfp);
    eigen_vectors = eigen_values_vectors.vectors;

    v1 = real(eigen_vectors[:, 1]);
    v2 = real(eigen_vectors[:,2]);
    v3 = imag(eigen_vectors[:, 3]);

    v1 = normalize(v1, norm);
    v2 = normalize(v2, norm);
    v3 = normalize(v3, norm);

    v1 = transpose(v1);
    v2 = transpose(v2);
    v3 = transpose(v3);

    A = transpose([v1; v2; v3]);

    return A;
end

function get_left_right_side(ϵ_box, counts_points, fp, A, dim)

    left_side = zeros(counts_points, dim);
    right_side = zeros(counts_points, dim);
    n = 1;
    α1 = 0.0;
    α_vector = [α1, 0.0, 0.0];
    while n <= counts_points

        α3 = -ϵ_box + 2 * ϵ_box * (n) / (counts_points);
        α_vector[3] = α3;

        α_vector[2] = - ϵ_box;
        point = fp + A * α_vector;
        left_side[n, :] = point;

        α_vector[2] = ϵ_box;
        point = fp + A * α_vector;
        right_side[n, :] = point;

        n+=1;
    end
    return left_side, right_side;
end

function get_up_down_side(ϵ_box, counts_points, fp, A, dim)

    down_side = zeros(counts_points, dim);
    up_side = zeros(counts_points, dim);
    n = 1;
    α1 = 0.0;
    α_vector = [α1, 0.0, 0.0];
    while n <= counts_points

        α2 = -ϵ_box + 2 * ϵ_box * (n) / (counts_points);
        α_vector[2] = α2;

        α_vector[3] = - ϵ_box;
        point = fp + A * α_vector;
        down_side[n, :] = point;

        α_vector[3] = ϵ_box;
        point = fp + A * α_vector;
        up_side[n, :] = point;

        n+=1;
    end
    return down_side, up_side;
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

function plot3_inter(data; indexs = [1, 2, 3], width = 500, height = 500, scales = [1, 1, 1], type = "lines")
        
    f = Figure(resolution = (width, height))
    ax = LScene(f[1, 1], show_axis = true)
    scale!(ax.scene, scales[1], scales[2], scales[3])
    if type == "lines"
        lines!(ax, data[:, indexs[1]], data[:, indexs[2]], data[:, indexs[3]], linewidth = 1.0)
    else
        scatter!(ax, data[:, indexs[1]], data[:, indexs[2]], data[:, indexs[3]], markersize = 5.0)
    end
    display(GLMakie.Screen(), f)
end

function plot3(data; indexs = [1, 2, 3], width = 500, height = 500)
        
    f = Figure(resolution = (width, height))
    ax = Axis3(f[1, 1])
    scatter!(ax, data[:, indexs[1]], data[:, indexs[2]], data[:, indexs[3]], linewidth = 1.0)
    display(GLMakie.Screen(), f)
end