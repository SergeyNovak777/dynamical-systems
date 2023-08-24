username = "Alex";
env = "bifurcation";
pathtorepo = "C:\\Users\\" *username *  "\\Desktop\\";
using Pkg;
Pkg.activate(pathtorepo * "dynamical-systems\\env\\" * env * "\\");


@inbounds function TM_bk(u, p)
    U(y, p) = p[8] + p[9] / ( 1.0 + exp( -50.0 * (y - p[7]) ) );
    σ(x, p) = 1.0 / ( 1.0 + exp( -20.0 * (x-p[6]) ) );
    g(E, x, y, p, U_) = log( 1.0 + exp( (p[5] * U_ * x * E + p[11]  ) / (p[1]) ) );
    
    U_ = U(u[3], p);
    du1 = (-u[1] + p[1] * g(u[1], u[2], u[3], p, U_) ) / p[2];
    du2 = (1.0 - u[2]) / p[3] - U_*u[2]*u[1];
    du3 = (-u[3])/p[4] + p[10] * σ(u[2], p);
    
    return [du1, du2, du3]
end

function main(){
    τ = 0.013;  τD = 0.07993;  τy = 3.3;  J = 3.07;  β = 0.300;
    xthr = 0.75; ythr = 0.4; α = 1.58; ΔU0 = 0.305;
    I0 = -1.6; U0 = 0.265;
    
    u0 = [8.638526524981895, 0.7320692774159767, 0.40718205935401675];
}


