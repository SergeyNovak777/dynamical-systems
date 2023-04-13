pathx= "I0_hom.mat"
pathy= "u0_hom.mat"

file = matopen(pathx)
I0_ = read(file, "I0_hom")
close(file)

file = matopen(pathy)
U0_ = read(file, "u0_hom")
close(file)

save("I0_hom.jld", "data",I0_); save("u0_hom.jld", "data",U0_)