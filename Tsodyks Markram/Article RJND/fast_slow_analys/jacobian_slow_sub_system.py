import sympy as sp
sp.init_printing()

x, y = sp.symbols('x, y', is_real = True)
taud, U0, deltaU0, ythr, tauy, beta, xthr, E = sp.symbols('tau_D, U0, Delta_U0, ythr, tau_y, beta, xtthr, E')

Uy = U0 + deltaU0 / ( 1 + sp.exp( -50 * ( y - ythr ) ) )
sigmay = 1 / ( 1 + sp.exp( -20 * ( x - xthr ) ) )
rhs_x =  ( 1 - x ) / taud - Uy * x * E;rhs_x
rhs_y = -y/tauy + beta * sigmay; rhs_y
rhs = sp.Matrix([rhs_x, rhs_y]); rhs

Jac = rhs.jacobian([x, y]); Jac