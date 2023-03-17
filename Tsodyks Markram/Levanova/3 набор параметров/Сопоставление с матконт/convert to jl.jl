pathx= "C:\\Users\\Alex\\Desktop\\x.mat"
pathy= "C:\\Users\\Alex\\Desktop\\y.mat"

file = matopen(pathx)
I0_ = read(file, "I0")
close(file)

file = matopen(pathy)
U0_ = read(file, "U0")
close(file)

save("I0.jld", "data",I0_); save("U0.jld", "data",U0_)