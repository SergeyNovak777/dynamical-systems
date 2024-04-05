I0_data = load("I0.jld")["data"]
U0_data = load("U0.jld")["data"]
I0_data = I0_data[:]
U0_data = U0_data[:];

function equat(x, c)
    c[1] .+ c[2]*x .+ c[3]*x.^2 .+ c[4]*x.^3 .+ c[5]*x.^4
end

xtest = equat(U0_data, f4.coeffs);

f = Figure(resolution = (900, 900))
axis = Axis(f[1, 1], xlabel = "I0", ylabel = "U0")
lines!(axis, xtest, U0_data, linewidth = 3.0, color = :lime)
scatter!(axis, I0_data, U0_data, markersize = 5.0, color = :black)
f