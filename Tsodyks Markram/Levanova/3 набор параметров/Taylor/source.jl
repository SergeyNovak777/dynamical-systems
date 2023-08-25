username = "Alex";
env = "integrate";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");

include(pathtorepo * "dynamical-systems\\system.jl");
using TaylorIntegration, GLMakie, JLD, StaticArrays, DynamicalSystems, LinearAlgebra;

function get_hom_curve()
    cd("C:\\Users\\Alex\\Desktop\\dynamical-systems\\Tsodyks Markram\\Levanova\\3 набор параметров\\Сопоставление с матконт\\файлы matlab");
    I0_hom = load("I0_hom_hom.jld")["data"];
    u0_hom = load("U0_hom_hom.jld")["data"];
    I0_hom = I0_hom[:];
    U0_hom = u0_hom[:];
    return I0_hom, U0_hom;
end

@inbounds function TM!(du, u, p, t)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) );
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) );
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) );
    
    U_ = U(u[3], p);
    du[1] = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) ) / p[2];
    du[2] = (1.0 - u[2]) / p[3] - U_*u[2]*u[1];
    du[3] = (-u[3])/p[4] + p[10] * σ(u[2], p);
end

function get_fp(p)

    intervalstart = -100; intervalend = 100;
    Ei = interval(intervalstart, intervalend);
    xi = interval(intervalstart, intervalend);
    yi = interval(intervalstart, intervalend);
    box = IntervalBox(Ei, xi, yi);

    ds = CoupledODEs(TM, [0.0, 0.0, 0.0], p);
    fp, ei, _ = fixedpoints(ds, box, jacob_TM_);

    return fp
end

function get_shift(p)

    fp = get_fp(p);
    index_fp = 1;
    index_vec = 1;
    ϵ = 1e-11;

    Jac = jacob_TM_(fp[index_fp], p, 0);
    eivecs = eigvecs(Jac);
    println(eivecs);
    shift =  fp[index_fp] + real(eivecs[:, index_vec])*ϵ;
    return shift

end
function main(order::Int64, ts::Int64 = 1, tf::Int64 = 50000)

    τ = 0.013;  τD = 0.07993;  τy = 3.3;  J = 3.07;  β = 0.300;
    xthr = 0.75; ythr = 0.4; α = 1.58; ΔU0 = 0.305;

    I0_hom, U0_hom = get_hom_curve()

    I0 = I0_hom[1]; U0 = U0_hom[1];

    p = [α, τ, τD, τy, J, xthr, ythr, U0, ΔU0, β, I0];
    tstart = 0.0; tfinish = -100.0;

    u0 = get_shift(p);
    u0 = Vector(u0);
    println(typeof(u0));
    
    t, sol = taylorinteg(TM!, u0, tstart, tfinish, order, 1e-20, p, maxsteps = 100000);

    indexx,indexy,indexz = 2, 3, 1

    GLMakie.activate!();
    f = Figure(resolution = (700, 700));
    ax = LScene(f[1, 1], show_axis = true);
    scale!(ax.scene, 50, 50, 1);
    lines!(ax, sol[ts:tf, indexx], sol[ts:tf, indexy], sol[ts:tf, indexz], linewidth = 1.5, color = :black);
    display(GLMakie.Screen(), f);

    return sol;

end


