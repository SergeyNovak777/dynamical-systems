function rulkov(u, p, n)

    function x_func(u, p)
        if u[1] <= 0
            return p[1] / ( 1 - u[1] ) + u[2]
        elseif u[1] > 0 && u[1] < p[1] + u[2]
            return p[1] + u[2]
        else
            return -1
        end
    end
        
    α, σ, μ = p
    x, y = u
    
    xn = x_func(u, p)
    yn = y + μ * ( -x - 1 + σ )
    
    return SVector{2}(xn, yn)
end

function rulkov_get_params()
    α = 5.6; σ = -0.25; μ = 0.001
    return [α, σ, μ]
end