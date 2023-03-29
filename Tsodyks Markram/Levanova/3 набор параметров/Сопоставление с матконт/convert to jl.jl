pathx= "I0doubleext.mat"
pathy= "U0doubleext.mat"

file = matopen(pathx)
I0_ = read(file, "I0doubleext")
close(file)

file = matopen(pathy)
U0_ = read(file, "U0doubleext")
close(file)

save("I0doubleext.jld", "data",I0_); save("U0doubleext.jld", "data",U0_)