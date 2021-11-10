# Dynamic FBA following the static optimization approach
# Mahadevan et al 2002
# Simple ecoli case study with 4 reactions and 3 external metabolites plus biomass

using Revise;
using MPCCLibrary;
using PyPlot;

## function for creating the dfba model
function dfba_base_model(A, dt ;kwargs...)
    ncp = get(kwargs, :ncp, 3);
    Adot = get(kwargs, :Adot, collocation_matrix(ncp));
    glc_id = 1;
    ac_id = 3;
    o2_id = 2;
    kLa = get(kwargs, :kLa, 7.5);       # h^-1
    vmax_O2 = get(kwargs, :vmax_O2, 15);    # mmol/gDW.h
    vmax_glc = get(kwargs, :vmax_glc, 10);   # mmol/gDW.h
    Km = get(kwargs, :Km, 0.015);      # mmol/L
    # vmax_dot = get(kwargs, :vmax_dot, [0.1; 0.3; 0.3; 0.1]);
    zglc_in = get(kwargs, :zglc_in, 5);   # mmol/L
    nfe = get(kwargs, :nfe, 1); 
    print_level = get(kwargs, :print_level, 0);

    nm, nr = size(A);
    model = create_model(linear_solver = "ma57", print_level = print_level);

    # model variables
    @variable(model, v[1:nfe, 1:nr] >= 0);
    @variable(model, z[1:nfe, 1:nm, 1:ncp+1] >= 0);   # dynamic
    @variable(model, x[1:nfe, 1:ncp+1] >= 0);         # dynamic
    @variable(model, aux[1:nfe, 1:nr]);               # auxiliary variable to avoid inequality with quadratic function
    @variable(model, D[1:nfe]);                 # dilution rate

    # model equality constraints
    @constraint(model, ode_zg[i in 1:nfe, j in [glc_id], k in 2:ncp+1], 
                sum(z[i,j,l]*Adot[l,k] for l in 1:ncp+1) - dt*(sum(A[j,r]*v[i,r] for r in 1:nr)*x[i,k] 
                + D[i]*(zglc_in - z[i,j,k])) == 0);
    @constraint(model, ode_za[i in 1:nfe, j in [ac_id], k in 2:ncp+1], 
                sum(z[i,j,l]*Adot[l,k] for l in 1:ncp+1) - dt*(sum(A[j,r]*v[i,r] for r in 1:nr)*x[i,k] 
                - D[i]*z[i,j,k]) == 0);
    @constraint(model, ode_zo[i in 1:nfe, j in [o2_id], k in 2:ncp+1], 
                sum(z[i,j,l]*Adot[l,k] for l in 1:ncp+1) - dt*(sum(A[j,r]*v[i,r] for r in 1:nr)*x[i,k] 
                + kLa*(0.21 - z[i, o2_id, k]) - D[i]*z[i,j,k]) == 0);
    @constraint(model, ode_x[i in 1:nfe, k in 2:ncp+1], 
                sum(x[i,l]*Adot[l,k] for l in 1:ncp+1) - dt*(sum(v[i,r] for r in 1:nr)*x[i,k] 
                - D[i]*x[i,k]) == 0);

    # continuity constraints
    if nfe > 1
        @constraint(model, cont_z[i in 2:nfe, j in 1:nm], z[i,j,1] - z[i-1,j,end] == 0);
        @constraint(model, cont_x[i in 2:nfe], x[i,1] - x[i-1,end] == 0);
    end

    # auxiliary constraint
    @constraint(model, aux_con[i in 1:nfe, r in 1:nr], aux[i,r] - v[i,r]*z[i,glc_id,end] == 0);

    # model inequality constraints
    @constraint(model, glc_lim_ub[i in 1:nfe], -sum(A[glc_id,r]*aux[i,r] for r in 1:nr) 
                - sum(A[glc_id,r]*v[i,r] for r in 1:nr)*Km + vmax_glc*z[i,glc_id,end] >= 0);
    @constraint(model, glc_lim_lb[i in 1:nfe], sum(A[glc_id,r]*aux[i,r] for r in 1:nr) 
                + sum(A[glc_id,r]*v[i,r] for r in 1:nr)*Km + vmax_glc*z[i,glc_id,end] >= 0);
    @constraint(model, o2_limit_ub[i in 1:nfe], sum(A[o2_id,r]*v[i,r] for r in 1:nr) + vmax_O2 >= 0);
    @constraint(model, ac_limit_ub[i in 1:nfe], -sum(A[ac_id,r]*v[i,r] for r in 1:nr) + 100.0 >= 0);

    return model;
end

function objective_model_static(model)
    v = model[:v];
    nfe, nr = size(v);
    @objective(model, Min, -sum(v[i,r] for i in 1:nfe for r in 1:nr));    # maximize growth
end

function update_model_static(model, z0, x0, v0)
    z = model[:z];
    x = model[:x];
    v = model[:v];
    fix.(z[1,:,1], z0, force = true);
    fix(x[1,1], x0, force = true);
    optimize!(model);
    return (value.(z[1,:,end]), value(x[1,end]), value.(v[1,:]));
end
