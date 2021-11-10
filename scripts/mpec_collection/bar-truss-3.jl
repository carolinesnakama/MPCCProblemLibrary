#
# MPEC Collection by Sven Leyffer
# Problem: bar-truss-3
#

# using Revise;
using MPCCLibrary;

model = create_model(print_level = 0);

# number of elements
nd = 2;    # structure dof
nm = 3;    # members
ny = 2;    # yield functs per member

# Parameters
L = [500; 400; 500];    # Element lengths
E = 20e3;               # Young's modulus
sigma = 50;             # Yield limit
F = [400; 600];
C = [ 0.6 0.8
      0.0 1.0
     -0.6 0.8];
N = hcat(ones(nm),-1*ones(nm));

# Variables
@variable(model, S[1:nm], start = 1);
@variable(model, r[1:nm, 1:ny], start = 1);
@variable(model, H[1:nm, 1:ny, 1:ny], start = 1);     # hardening parameters in tension and compression
@variable(model, Q[1:nm]);
@variable(model, -4 <= u[1:nd] <= 4);                 # deflection
@variable(model, a[1:nm] >= 0, start = 1);            # bar areas
@variable(model, z[1:nm, 1:ny] >= 0);
@variable(model, w[1:nm, 1:ny] >= 0);                 # yield function

# Constraints
@constraint(model, tech[i in 1:nm], a[i] == a[1]);
@constraint(model, stiff[i in 1:nm], S[i] == E*a[i]/L[i]);
@constraint(model, limit[i in 1:nm, j in 1:ny], r[i,j] == sigma*a[i]);
@constraint(model, hard[i in 1:nm, j in 1:ny], H[i,j,j] == 0.125*S[i]);
@constraint(model, compat[i in 1:nm], Q[i] == S[i]*sum(C[i,k] * u[k] for k in 1:nd)
                                        - S[i]*sum( N[i,j] * z[i,j] for j in 1:ny));
@constraint(model, equil[k in 1:nd], F[k] == sum(C[i, k]*Q[i] for i in 1:nm));
@constraint(model, yield[i in 1:nm, j in 1:ny], w[i,j] == -N[i,j]*Q[i] 
                                                    + sum(H[i,j,jj]*z[i,jj] for jj in 1:ny) + r[i,j]);                                        

# define original objective function as an expression
@expression(model, obj, sum(L[i]*a[i] for i in 1:nm));

# H is a diagonal matrix - setting off-diagonal values to zero 
fix.(H[:,1,2], 0.0);
fix.(H[:,2,1], 0.0);

# Complentarity pairs
comp = [(w, z)];

# solving the MPCC problem
ans = solve_mpcc!(model, comp, false, :obj);