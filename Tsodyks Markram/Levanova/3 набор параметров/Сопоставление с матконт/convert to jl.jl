using Pkg
Pkg.activate("C:\\Users\\Alex\\Desktop\\dynamical-systems\\env\\bifurcation\\")
cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab")
using MAT, JLD
pathx= "vars.mat"
#pathy= "U0_hom_bt.mat"

file = matopen(pathx)
I0_ = read(file, "vars")
close(file)

#file = matopen(pathy)
#U0_ = read(file, "U0_hom_bt")
#close(file)

save("vars.jld", "data",I0_); #save("U0_hom_bt.jld", "data",U0_)