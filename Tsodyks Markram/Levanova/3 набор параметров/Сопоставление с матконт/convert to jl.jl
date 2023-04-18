using Pkg
Pkg.activate("C:\\Users\\Alex\\Desktop\\dynamical-systems\\env\\bifurcation\\")
cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab")
using MAT, JLD
pathx= "I0_hom.mat"
pathy= "u0_hom.mat"

file = matopen(pathx)
I0_ = read(file, "I0_hom")
close(file)

file = matopen(pathy)
U0_ = read(file, "u0_hom")
close(file)

save("I0_hom.jld", "data",I0_); save("u0_hom.jld", "data",U0_)