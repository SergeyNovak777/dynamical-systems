pathx= "I0_double_detailed_1.mat"
pathy= "U0_double_detailed_1.mat"

file = matopen(pathx)
I0_ = read(file, "I0_double_")
close(file)

file = matopen(pathy)
U0_ = read(file, "U0_double_")
close(file)

save("I0_double_detailed_1.jld", "data",I0_); save("U0_double_detailed_1.jld", "data",U0_)