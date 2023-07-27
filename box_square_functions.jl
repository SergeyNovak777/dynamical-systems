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

function get_norm_αs(u, fp, A)
    x = Vector(u)
    b = x - fp
    b = Vector(b)
    linprob = LinearProblem(A, b)
    linsolve = solve(linprob)
    return norm(linsolve.u, Inf), linsolve
end

function get_fixed_point(system, jac_system, p, u0)
    ds = CoupledODEs(system, u0, p)
    fp, _, _ = fixedpoints(ds, box, jac_system)
    return fp
end

function get_matrix(fixedpoint, p, jac_system, ϵ_box, norm = 2)
    jac_fp = jac_system(fixedpoint, p, 0.0)
    
    eigen_val_vec = eigen(jac_fp)
    eigen_vectors = eigen_val_vec.vectors

    v1 = real(eigen_vectors[:, 1])
    v2 = real(eigen_vectors[:, 2])
    v3 = imag(eigen_vectors[:, 3])

    v1 = normalize(v1, norm)
    v1t = transpose(v1)

    v2 = normalize(v2, norm)
    v2t = transpose(v2)

    v3 = normalize(v3, norm)
    v3t = transpose(v3)

    A = transpose([v1t; v2t; v3t])
    A = Matrix(A)
end

function get_arrays_dots(number_points_on_side, dim)
        dots_right_side = zeros(number_points_on_side, dim)
        array_α_vec_right = zeros(number_points_on_side, dim)
    
        dots_left_side = zeros(number_points_on_side, dim)
        array_α_vec_left = zeros(number_points_on_side, dim)
    
        dots_up_side = zeros(number_points_on_side, dim)
        array_α_vec_up = zeros(number_points_on_side, dim)
    
        dots_down_side = zeros(number_points_on_side, dim)
        array_α_vec_down = zeros(number_points_on_side, dim)
    
        return dots_right_side, array_α_vec_right,
                dots_left_side, array_α_vec_left,
                dots_up_side, array_α_vec_up,
                dots_down_side, array_α_vec_down
end

function fill_side_square(ϵ_box, fixedpoint, number_points_on_side, A,
    dots_right_side, array_α_vec_right, dots_left_side, array_α_vec_left,
    dots_up_side, array_α_vec_up, dots_down_side, array_α_vec_down)
    
    # case α1 = 0; α2 = ϵ_box; α3 = [-ϵ_box; ϵ_box]
    n = 1
    N = number_points_on_side
    while n <= number_points_on_side

        α1 = 0.0
        α2 = ϵ_box
        α3 = - ϵ_box + 2 * ϵ_box * n / N
        α_vec = [α1, α2, α3]
        point = fixedpoint + A * α_vec
        
        dots_right_side[n, :] = point
        _, α_vec = get_norm_αs( dots_right_side[n, :], fixedpoint, A )
        array_α_vec_right[n, :] = α_vec
        n+=1
    end

    # case α1 = 0; α2 = -ϵ_box; α3 = [-ϵ_box; ϵ_box]
    n = 1
    N = number_points_on_side
    while n <= number_points_on_side

        α1 = 0.0
        α2 = -ϵ_box
        α3 = - ϵ_box + 2 * ϵ_box * n / N
        α_vec = [α1, α2, α3]
        point = fixedpoint + A * α_vec
        
        dots_left_side[n, :] = point
        _, α_vec = get_norm_αs( dots_left_side[n, :], fixedpoint, A )
        array_α_vec_left[n, :] = α_vec
        n+=1
    end

    # case α1 = 0; α2 = [-ϵ_box; ϵ_box]; α3 = ϵ_box
    n = 1
    N = number_points_on_side
    while n <= number_points_on_side

        α1 = 0.0
        α2 = - ϵ_box + 2 * ϵ_box * n / N
        α3 = ϵ_box
        α_vec = [α1, α2, α3]
        point = fixedpoint + A * α_vec
        
        dots_up_side[n, :] = point
        _, α_vec = get_norm_αs( dots_up_side[n, :], fixedpoint, A )
        array_α_vec_up[n, :] = α_vec
        n+=1
    end

    # case α1 = 0; α2 = [-ϵ_box; ϵ_box]; α3 = -ϵ_box
    n = 1
    N = number_points_on_side
    while n <= number_points_on_side

        α1 = 0.0
        α2 = - ϵ_box + 2 * ϵ_box * n / N
        α3 = -ϵ_box
        α_vec = [α1, α2, α3]
        point = fixedpoint + A * α_vec
        
        dots_down_side[n, :] = point
        _, α_vec = get_norm_αs( dots_down_side[n, :], fixedpoint, A )
        array_α_vec_down[n, :] = α_vec
        n+=1
    end
end

function trajectory_from_side(p, A, dots_from_side, number_points, fixedpoint,  cb, data)
    data = [dots_u0, check_events, time_events, dots_on_event, αs, norms]
    for index in range(1, number_points, step = 1)
    
        u0 = SA[dots_from_side[index, 1], dots_from_side[index, 2], dots_from_side[index, 3]]
        prob = ODEProblem(TM, u0, tspan, p)
        sol = solve(prob, alg = Vern9(), abstol = 1e-14, reltol = 1e-14, callback = cb)
        norm_, linsolve = get_norm_αs(sol[end], fixedpoint, A)
        if sol.retcode == ReturnCode.Terminated
            dots_u0[index, 1] = sol[end][1]
            dots_u0[index, 2] = sol[end][2]
            dots_u0[index, 3] = sol[end][3]
            check_events[index] =  true
            time_events[index] = sol.t[end]
            αs[index, 1] = linsolve[1]
            αs[index, 2] = linsolve[2]
            αs[index, 3] = linsolve[3]
            norms[index] = norm_
        else
            check_events[index] =  false
        end
    end

end