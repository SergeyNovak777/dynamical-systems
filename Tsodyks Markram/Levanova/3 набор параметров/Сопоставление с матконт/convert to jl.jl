using Pkg
Pkg.activate("C:\\Users\\Alex\\Desktop\\dynamical-systems\\env\\bifurcation\\")
cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab")
using MAT, JLD
pathx= "I0_hom_stump.mat"
pathy= "U0_hom_stump.mat"

file = matopen(pathx)
I0_ = read(file, "I0_hom_stump")
close(file)

file = matopen(pathy)
U0_ = read(file, "U0_hom_stump")
close(file)

save("I0_hom_stump.jld", "data",I0_); save("U0_hom_stump.jld", "data",U0_)