versioninfo()

Threads.nthreads()

Threads.threadid()

for i in 1:Threads.nthreads()
    println("i: $(i) \t Thread ID: $(Threads.threadid())")
end

Threads.@threads for i in 1:Threads.nthreads()
    println("i: $(i) \t Thread ID: $(Threads.threadid())")
end