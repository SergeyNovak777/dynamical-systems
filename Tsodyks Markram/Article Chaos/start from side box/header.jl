function get_fixed_point(sys, jac_sys, params, dim, box)
    ds = CoupledODEs(sys, zeros(dim), params)
    fp, _, _ = fixedpoints(ds, box, jac_sys)
end

function get_matrix(jac_sys, params, fixedpoint, norm_)

    Ju0 = jac_sys(fixedpoint, params, 0.0);
    eigen_vec_val = eigen(Ju0);
    eigen_vec = eigen_vec_val.vectors;

    v1 = eigen_vec[:, 1];
    v2 = eigen_vec[:, 2];
    v3 = eigen_vec[:, 3];

    v1 = real(v1);
    v2 = real(v2);
    v3 = imag(v3);

    v1 = normalize(v1, norm_);
    v2 = normalize(v2, norm_);
    v3 = normalize(v3, norm_);

    v1 = transpose(v1);
    v2 = transpose(v2);
    v3 = transpose(v3);

    A = transpose([v1; v2; v3]);
    A = Matrix(A);
    return A;
end

function get_event(fixedpoint, A, ϵ_box)
    function condition(u, p ,t)
        u = Vector(u);
        b = u - fixedpoint;
        b = Vector(b);
        linprob = LinearProblem(A, b);
        linsolve = solve(linprob);
        return norm(linsolve.u, Inf) - ϵ_box;
    end
    return condition;
end

affect!(integrator) = terminate!(integrator)