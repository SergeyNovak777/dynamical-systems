using Pkg
Pkg.activate("C:\\Users\\Alex\\Desktop\\dynamical-systems\\env\\bifurcation\\")
cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab")
using MAT, JLD
pathx= "x_to_top_righter.mat"
#pathy= "U0_hom_bt.mat"

file = matopen(pathx)
output = read(file, "x")
close(file)

#file = matopen(pathy)
#U0_ = read(file, "U0_hom_bt")
#close(file)

save("x_to_top_righter.jld", "data",output); #save("U0_hom_bt.jld", "data",U0_)