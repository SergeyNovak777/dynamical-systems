if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\implicit\\")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/implicit/")
end

using Implicit3DPlotting

f(x) = x[3] - 1.58 * log( 1 + exp( ( 3.07 * ( 0.265 + 0.305 / ( 1 + exp( -50 * ( x[2]  -0.4 ) ) ) ) * x[1] * x[3] -1.6  ) / 1.58 ) )

scene = plot_implicit_surface(f; transparency=false, xlims=(0.5,1.0), ylims=(0.4,1.0), zlims=(0,25))

