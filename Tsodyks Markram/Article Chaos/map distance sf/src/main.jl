include("/home/sergey/work/repo/dynamical-systems/Tsodyks Markram/Levanova/3 набор параметров/map distance sf/src/header.jl");

function main()
    
    time_integrate = 5000;
    tstep = 0.001;
    integ_set = (alg = RK4(), adaptive = false, dt = tstep);

    τ = 0.013;  τD = 0.07993;  τy = 3.3;  J = 3.07;  β = 0.300
    xthr = 0.75; ythr = 0.4
    α = 1.58; ΔU0 = 0.305;

    p = [α, τ, τD, τy, J, xthr, ythr, 0.0, ΔU0, β, 0.0]

    len = 400;
    I0range = range( -1.7, -1.73,  length = 400 )
    U0range = range( 0.268,  0.265, length = 400 )

    global dis_s = zeros(length(I0range), length(U0range));
    
    path_to_file = "/home/sergey/work/repo/dynamical-systems/Tsodyks Markram/Levanova/3 набор параметров/Сопоставление с матконт/Карты спектров/"
    u0s = load(path_to_file*"u0s_400_400_tauD_article2_zoom_for_curve.jld")["data"];

    lenI0 = length(I0range)
    lenU0 = length(U0range)

    map_dim = " $(lenI0)x$(lenU0) ";
    name = " tauD_article2_map_distance";
    format = ".jld2";
    namefile_dis = "LSE" * map_dim * name * format;
    # Индексы фиксируемого и управляющего параметра
    # для предварительной протяжки
    index_fix = 11;
    var_fix = I0range[1];
    index_control = 8;

    # Индексы управляющих параметров
    index_p2 = 11; #I0
    index_p1 = 8; #U0
    
    for (idx_U0, U0_) in enumerate(U0range)
        for (idx_I0, I0_) in enumerate(I0range)
            
            u0_lc = u0s[idx_I0, idx_U0, :]

            ds = init_ds(p, index_p1, index_p2, U0_, I0_, u0_lc, integ_set)
            dis = calculate_distance(ds, time_integrate)
            
            #output(idx_I0, idx_U0, I0_, U0_, u0_lc, dis)
            #separate()

            save_output(idx_I0, idx_U0, dis)
        end
        save_tofile(namefile_dis)
    end
end