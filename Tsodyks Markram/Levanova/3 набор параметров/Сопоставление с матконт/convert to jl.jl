using Pkg
Pkg.activate("C:\\Users\\Alex\\Desktop\\dynamical-systems\\env\\bifurcation\\")
cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab")
using MAT, JLD
pathx= "I0_Fold.mat"
pathy= "U0_Fold.mat"

file = matopen(pathx)
I0_ = read(file, "I0_Fold")
close(file)

file = matopen(pathy)
U0_ = read(file, "U0_Fold")
close(file)

save("I0_Fold.jld", "data",I0_); save("U0_Fold.jld", "data",U0_)