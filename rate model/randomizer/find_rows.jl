
if Sys.iswindows()
    username = "Alex"
    pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\"
    using Pkg
    Pkg.activate(pathtorepo * "dynamical-systems\\env\\integrate\\")
else
    username = "sergey"
    pathtorepo = "/home/" *username *"/work/repo/dynamical-systems"
    using Pkg
    Pkg.activate(pathtorepo * "/env/integrate/")
end

using DataFrames, CSV

df = DataFrame(CSV.File("/home/sergey/work/repo/dynamical-systems/rate model/randomizer/50000_u0s.csv"))
df_filter = filter(r -> r.LE1 >= 0, df)

u0s = df_filter[:,[:"sE_start",:"sI_start", :"rE_start", :"rI_start", :"Y_start"]]