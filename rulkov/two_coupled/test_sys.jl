if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
end

using StaticArrays, DynamicalSystems

function two_coupled_rulkov_first_iteration(u, p)

    function right_part_x(x, y, z)
        if x <= 0.0
            return α/(1.0 - x) + y
        elseif (0.0 < x < (alpha + y)) && (z <=0.0)
            return α  +y
        elseif (x >= α + y) || (z > 0)
            return -1.0
        else
            return -1.0;
        end
    end

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, I21prev, x2, y2, z2, I12prev = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, k, γ_1, γ_2, g1, g2 = p

    I21 = γ_2 * I21prev + g2 * ( x_rp - x1 ) * xi(x2)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21, z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21)
    z1n = x1

    I12 = γ_1 * I12prev + g1 * ( x_rp - x2 ) * xi(x1)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12, z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12)
    z2n = x2

    r

function two_coupled_rulkov(u, p, t)

    function right_part_x(x, y, z)
        if x <= 0.0
            return α/(1.0 - x) + y
        elseif (0.0 < x < (alpha + y)) && (z <=0.0)
            return α  +y
        elseif (x >= α + y) || (z > 0)
            return -1.0
        else
            return -1.0;
        end
    end

    function right_part_y(x, y, sum_inhibitory)
        return y + μ * ( -x - 1.0 + σ + sum_inhibitory )
    end

    function xi(x)
        if x > x_th
            return 1.0;
        elseif x <=x_th
            return 0.0;
        end
    end

    x1, y1, z1, I21prev, x2, y2, z2, I12prev = u
    α, σ, μ, β_syn, σ_syn, x_rp, x_th, k, γ_1, γ_2, g1, g2 = p

    I21 = γ_2 * I21prev + g2 * ( x_rp - x1 ) * xi(x2)
    x1n = right_part_x(x1, y1 + (β_syn/k) * I21, z1)
    y1n = right_part_y(x1, y1, (σ_syn/k) * I21)
    z1n = x1

    I12 = γ_1 * I12prev + g1 * ( x_rp - x2 ) * xi(x1)
    x2n = right_part_x(x2, y2 + (β_syn/k) * I12, z2)
    y2n = right_part_y(x2, y2, (σ_syn/k) * I12)
    z2n = x2

    return [x1n, y1n, z1n, I21, x2n, y2n, z2n, I12]
end

function get_params_two_coupled_rulkov()
    α = 4.6; σ =-0.1; μ = 0.001;
    β_syn = 0.0001; σ_syn = 1.0; k = 1.0;
    x_rp = -1.5; x_th = -0.8;
    γ_1 = 0.0; γ_2 = 0.0;
    g1 = 0.8; g2 = 0.2;
    return [α, σ, μ, β_syn, σ_syn, x_rp, x_th, k, γ_1, γ_2, g1, g2]
end

