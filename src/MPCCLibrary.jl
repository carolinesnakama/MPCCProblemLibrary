module MPCCLibrary

using Ipopt, LinearAlgebra, Printf;
using Combinatorics, Polynomials;

using Reexport;
    @reexport using JuMP

export 
    
        create_model, solve_mpcc!, 

        collocation_matrix, collocation_mat;

include("mpcc.jl");

include("collocation.jl");

end