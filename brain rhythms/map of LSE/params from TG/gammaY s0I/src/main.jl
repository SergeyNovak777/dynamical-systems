include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\brain rhythms\\map of LSE\\params from TG\\gammaY s0I\\src\\header.jl");

#include("C:\\Users\\Alex\\Desktop\\dynamical-systems\\brain rhythms\\map of LSE\\params from TG\\gammaY s0I\\src\\main.jl")
function main()

    time_LSE = 1000;
    time_attract = 1000;
    tstep = 0.001;
    integ_set = (alg = RK4(), adaptive = false, dt = tstep);

    τsE = 3.0; γE = 4.0; s0E = 0.15;
    τsI = 10.0; γI = 8.0; s0I = 0.1;
    τrE = 2.0; kE = 5.0; IE = 0.9; wEE = 3.5; wIE = 5.0; θE = 0.2;
    τrI = 6.0; kI = 5.0; II = 0.0; wEI = 5.0; wII = 3.0; θI = 0.4;
    τY = 0.01;  βY = 0.01;
    ythr = 0.5; sEthr = 0.5; kY = 0.01
    γY = 0.0;

    p = [τsE, γE, s0E, τsI, γI, s0I, τrE, kE, IE, wEE, wIE, θE, τrI, kI, II,
        wEI, wII, θI, τY, βY, γY, ythr, sEthr, kY];
    # τsE - 1, γE - 2, s0E - 3, τsI - 4, γI - 5, s0I - 6, τrE - 7, kE - 8, IE - 9,
    # wEE - 10, wIE - 11, θE - 12, τrI - 13, kI - 14, II - 15,
    # wEI - 16, wII - 17, θI - 18, τY - 19, βY - 20, γY - 21, ythr - 22, sEthr - 23, kY - 24

    dim = 5
    u0 = zeros(dim);

    len = 100;
    p1name = "γY";
    p2name = "s0I";
    p1_range = range( 0.0, 10.0, length = len );
    p2_range = range(0.0, 0.5, length = len);

    global Λs = zeros(len, len, dim);
    global u0s = zeros(len, len, dim);

    map_dim = " $(len)x$(len) ";
    name = " $(p1name) $(p2name) rate_model params from tg RK4";
    format = ".jld2";
    namefile_LSE = "LSE" * map_dim * name * format;
    namefile_u0s = "u0s" * map_dim * name * format;

    # Индексы фиксируемого и управляющего параметра
    
    # Индексы управляющих параметров
    index_p1= 21;
    index_p2 = 6;

    # для предварительной протяжки
    index_fix = index_p1;
    var_fix = p1_range[1];
    index_control = index_p2;
    
    for (p2_loc_index, p2_loc) in enumerate(p2_range)
    
        if p2_loc_index == 1
            global u0_lc = u0
        end
        
        #output(p2name, p2_loc_index,p2_loc, u0_lc)
        
        ds = init_ds_(rate_model, p, index_control, p2_loc,
        index_fix, var_fix, u0_lc, integ_set)
        #println("P $(p)");flush(stdout)
        u0_lc = goto_attractor(ds, time_attract, integ_set)

        ds = init_ds_(rate_model, p, index_control, p2_loc,
        index_fix, var_fix, u0_lc, integ_set)
        
        ΛΛ = spectrum(ds, time_LSE)
        
        #output_end_iter(ΛΛ, u0_lc)
        
        save_output(p2_loc_index, ΛΛ, u0_lc)
        save_tofile(namefile_LSE, namefile_u0s)
        #separate()
    end
   
    #println("  "); flush(stdout)
    for (p2_loc_index, p2_loc) in enumerate(p2_range)
        for (p1_loc_index, p1_loc) in enumerate(p1_range)
            
            if p1_loc_index == 1
                continue
            end
            
            u0_lc = u0s[p1_loc_index - 1, p2_loc_index, :]

            #output(p1name, p2name, p1_loc_index, p2_loc_index, p1_loc, p2_loc, u0_lc)
            
            ds = init_ds(rate_model, p, index_p2, index_p1, p2_loc, p1_loc, u0_lc, integ_set)
            #println("P $(p)");flush(stdout)
            u0_lc = goto_attractor(ds, time_attract, integ_set)
            ds = init_ds(rate_model, p, index_p2, index_p1, p2_loc, p1_loc, u0_lc, integ_set)
            ΛΛ = spectrum(ds, time_LSE)
            
            #output_end_iter(ΛΛ, u0_lc)
            
            save_output(p1_loc_index, p2_loc_index, ΛΛ, u0_lc)
            
            #separate()
        end
        save_tofile(namefile_LSE, namefile_u0s)
    end
end