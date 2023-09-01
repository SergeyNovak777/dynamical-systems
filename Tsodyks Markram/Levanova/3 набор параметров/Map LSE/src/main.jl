include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Map LSE\\src\\header.jl");

function main()

    time_LSE = 500;
    time_attract = 600;
    tstep = 0.001;
    integ_set = (alg = RK4(), adaptive = false, dt = tstep);

    τ = 0.013;  τD = 0.07993;  τy = 3.3;  J = 3.07;  β = 0.300
    xthr = 0.75; ythr = 0.4
    α = 1.58; ΔU0 = 0.305;

    p = [α, τ, τD, τy, J, xthr, ythr, 0.0, ΔU0, β, 0.0]
    u0 = [10.870367054955267, 0.6670395183801261, 0.44434050730193664];

    len = 400;
    I0range = range( -0.5, -2.0,  length = len );
    U0range = range( 0.15,  0.55, length = len );

    global Λs = zeros(length(I0range), length(U0range), 3);
    global u0s = zeros(length(I0range), length(U0range), 3);

    lenI0 = length(I0range)
    lenU0 = length(U0range)

    map_dim = " $(lenI0)x$(lenU0) ";
    name = " tau_D article2 first map extended";
    format = ".jld2";
    namefile_LSE = "LSE" * map_dim * name * format;
    namefile_u0s = "u0s" * map_dim * name * format;

    # Индексы фиксируемого и управляющего параметра
    # для предварительной протяжки
    index_fix = 11;
    var_fix = I0range[1];
    index_control = 8;

    # Индексы управляющих параметров
    index_p2 = 11; #I0
    index_p1 = 8; #U0
    

    for (idx_U0, U0_) in enumerate(U0range)
    
        if idx_U0 == 1
            global u0_lc = u0
        end
        
        #output(idx_U0,U0_, u0_lc)
        
        ds = init_ds_(p, index_control, U0_,
        index_fix, var_fix, u0_lc, integ_set)

        u0_lc = goto_attractor(ds, time_attract, integ_set)

        ds = init_ds_(p, index_control, U0_,
        index_fix, var_fix, u0_lc, integ_set)
        
        ΛΛ = spectrum(ds, time_LSE)
        
        #output_end_iter(ΛΛ, u0_lc)
        
        save_output(idx_U0, ΛΛ, u0_lc)
        save_tofile(namefile_LSE, namefile_u0s)
        #separate()
    end
   

    for (idx_U0, U0_) in enumerate(U0range)
        for (idx_I0, I0_) in enumerate(I0range)
            
            if idx_I0 == 1 # Если использовать while условие не нужно!
                continue
            end
            
            
            u0_lc = u0s[idx_I0 - 1, idx_U0, :]
            
            #output(idx_I0, idx_U0, I0_, U0_, u0_lc)
            
            ds = init_ds(p, index_p1, index_p2, U0_, I0_, u0_lc, integ_set)
            u0_lc = goto_attractor(ds, time_attract, integ_set)
            ds = init_ds(p, index_p1, index_p2, U0_, I0_, u0_lc, integ_set)
            ΛΛ = spectrum(ds, time_LSE)
            
            #output_end_iter(ΛΛ, u0_lc)
            
            save_output(idx_I0, idx_U0, ΛΛ, u0_lc)
            
            #separate()
        end
        save_tofile(namefile_LSE, namefile_u0s)
    end
end

"""
α -- 1, τ -- 2, τD -- 3, τy -- 4, J -- 5,
xthr -- 6, ythr -- 7, U0 -- 8, ΔU0 -- 9, β -- 10
I0 -- 11
"""