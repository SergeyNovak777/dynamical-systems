U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) )
σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) )
function TM!(du, u, p, t)
    du[1] = (-u[1] + p[1] * log( 1.0 + exp( (p[5] * U(u[3], p) * u[2] * u[1] + p[11]  ) / (p[1]) ) ) ) / p[2]
    du[2] = (1.0 - u[2])/p[3] - U(u[3], p)*u[2]*u[1]
    du[3] = (-u[3])/p[4] + p[10] * σ(u[2], p)
    return du
end
const τ = 0.013;  const τD = 0.080;  const τy = 3.3;  const J = 3.07;  const β = 0.300
const xthr = 0.75; const ythr = 0.4
const α = 1.58;  const U0 = 0.30;  const ΔU0 = 0.305;
I0 = -1.0

TM(u, p) = TM!(similar(u), u, p, 0)
p = (α, τ, τD, τy, J, xthr, ythr, U0, ΔU0, β, I0)
u0 = [9.25834,  0.732799,  0.410686]

prob = BifurcationProblem(TM, u0, p, (@lens _[11]))

opts_new = NewtonPar(maxIter = 3, tol = 1e-10)
opts_con = ContinuationPar(pMin = -3.5, pMax = 0.0,
                            ds = 0.00001, dsmin = 1e-10, dsmax = 0.001,
                            nev = 3, detectBifurcation = 3, newtonOptions  = opts_new,
                            maxSteps  = 5000)
br = continuation(prob, PALC(), opts_con;
                bothside = true,
                detectFold = true)



opts_new_ = NewtonPar(maxIter = 10, tol = 1e-6)
opts_con_ = ContinuationPar(pMin = -3.5, pMax = 0.0,
                            ds = 0.01, dsmin = 1e-5, dsmax = 0.1,
                            nev = 3, detectBifurcation = 3, newtonOptions  = opts_new_,
                            maxSteps  = 400)
Mt = 40
br_po = continuation( br, 4, opts_con_, PeriodicOrbitOCollProblem(Mt, 5))
#scene = plot(br_po, plotfold=false, markersize=2, legend=:topright)