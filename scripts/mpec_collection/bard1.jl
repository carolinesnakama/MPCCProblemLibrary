#
# MacMPEC Collection by Sven Leyffer
# Problem: bard1
#

# using Revise;
using MPCCLibrary;

model = create_model(print_level = 0);
@variable(model, x >= 0, start = 0);
@variable(model, y >= 0, start = 0);
@variable(model, l[1:3] >= 0, start = 0);   # lagrangian multipliers for of the inner optimization problem 

# constraint corresponding to condition ∇L = 0
@constraint(model, kkt, 2*(y - 1) - 1.5*x + l[1] - 0.5*l[2] + l[3] == 0);

# auxiliary variable for creating an array of expressions
tmp = Array{Any,1}(undef, 3);
tmp[1] = @expression(model, 3*x - y - 3);
tmp[2] = @expression(model, - x + 0.5*y + 4);
tmp[3] = @expression(model, - x - y + 7);

# constraints of the inner problem are defined as expressions for defining the complementarity term
@expression(model, g[i = 1:3], tmp[i]);
# inequality constraints of the inner problem
@constraint(model, posg[i = 1:3], g[i] >= 0);

# define original objective function as an expression
@expression(model, obj, (x - 5)^2 + (2*y + 1)^2);

ans = solve_mpcc!(model, [(g, l)], false, :obj);