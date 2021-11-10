#
# Functions for orthogonal collocation 
#
# Caroline Nakama
# October 2020
#
# Last updated: 23.11.2020
#

# Radau and Legendre roots for degrees 1 to 5
proots = Dict();
proots["Radau"] = [[1.],
                  [0.333333; 1.],
                  [0.155051; 0.644049; 1.],
                  [0.088588; 0.409467; 0.787659; 1.],
                  [0.057104; 0.276843; 0.583590; 0.860240; 1.]];
proots["Legendre"] = [[.5],
                     [0.211325; 0.788675],
                     [0.112702; 0.5; 0.887298],
                     [0.069432; 0.330009; 0.669991; 0.930568],
                     [0.046910; 0.230765; 0.5; 0.769235; 0.953090]];


"""
    collocation_mat(ncp = 3)

Returns the collocation matrix using `ncp` Radau roots as collocation points.

Note: DEPRECATED
"""
function collocation_mat(ncp = 3)                          # TODO: Implement it in the standard way
    if ncp != 3
        error("This function has only been implemented for 3 collocation points so far.")
    end
    radau = [0.15505, 0.64495, 1.00000];                       # TODO: More points
    M1 = zeros(ncp, ncp);
    M2 = zeros(ncp, ncp);
    for i in 1:ncp
        for j in 1:ncp
            M1[i,j] = radau[i]^j;
            M2[i,j] = j*radau[i]^(j-1);
        end
    end
    return M1*M2^-1;
end

"""
    collocation_matrix(ncp = 3, method = "Radau")

Function for generating a collocation matrix using `ncp` collocation points with `method` roots.

Returns the collocation matrix.

Possible values are 1 to 5 for `ncp`, and "Radau" and "Legendre" for `method`.

Reference: Lorenz T. Biegler. _Nonlinear programming: concepts, algorithms, and applications to chemical processes_. Vol. 10. Siam, 2010 (Section 10.2).
"""
function collocation_matrix(ncp = 3, method = "Radau")
    if ncp > 5
        error("This function can only take these values $(length.(proots["Radau"])) for ncp")
    end
    if method ∉ keys(proots)
        error("This method is not available. Possible values are $(keys(proots))")
    end
    t = vcat([0.], proots[method][ncp]);
    l = Array{Any,1}(undef, ncp+1);
    dl = copy(l);
    adot = zeros(ncp+1, ncp+1);
    comb = collect(combinations(Array(1:ncp+1), ncp));
    for i in 1:ncp+1
        k = setdiff(Array(1:ncp+1), comb[i])[1];
        l[k] = fromroots(t[comb[i]])/prod((t[k] - t[comb[i]][j]) for j in 1:ncp);
        dl[k] = derivative(l[k]);
        for j in 1:ncp+1
            adot[k,j] = dl[k](t[j]);
        end
    end
    return adot;
end