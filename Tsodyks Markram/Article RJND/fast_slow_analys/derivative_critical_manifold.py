import sympy as sp
sp.init_printing()

E, x, y = sp.symbols('E, x, y', is_real = True)
alpha, J, U0, deltaU0, ythr, I0 = sp.symbols('alpha, J, U_0, Delta_U_0, ythr, I0')

#alpha = 1.58; J = 3.07; U0 = 0.265; deltaU0 = 0.305; ythr = 0.4; I0 = -1.6;

Uy = U0 + deltaU0 / ( 1 + sp.exp( -50 * ( y - ythr ) ) )

f = -E + alpha * sp.log( 1 + sp.exp( ( J*Uy * x * E + I0 ) / alpha ) ); f

dfdE = sp.diff(f, E); dfdE
#dfdE.subs({x:0.441176, y:0.485294, E:0.642736})

sp.solve([f, dfdE],[ E,x,y])