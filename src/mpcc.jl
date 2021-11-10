#
# Implementation of Algorithm I from Leyffer et al 2005 
# Interior Methods for Mathematical Programs with Complementarity Constraints
#
# Caroline S. M. Nakama
# April 2020
#
# Last update: 28.10.2020
#

"""
    create_model(args...)

Returns a JuMP model using Ipopt as the optimization solver. `args` are attributes taken by the Ipopt optimizer 
(e.g. `print_level = 0`, `linear_solver = "ma57"`, etc. A list of all possible options can be found at the 
[Ipopt Documentation](https://coin-or.github.io/Ipopt/OPTIONS.html)).
"""
function create_model(;kwargs...)
    model = Model(Ipopt.Optimizer);
    if !isempty(kwargs)
        for (key, value) in kwargs
            set_optimizer_attributes(model, String(key) => value);
        end
    end
    return model;
end

"""
    solve_mpcc!(model, comp, plot_data = false, objective = :obj, args...)

Solves an MPCC problem defined in JuMP using the Ipopt solver.

`model` is a JuMP model initialized Ipopt Optimizer and containing the MPCC problem;
`comp` is an array of pairs with the complementarity variables;
`plot_data` is a boolean that defines whether the function returns data with the algorithm progress;
`objective` is a symbol corresponding to the name of the JuMP expression that defines the objective function.
`args` are parameters used by the MPCC algorithm.

Returns 0 (if `plot_data` is false) or the progress data from the algorithm (if `plot_data` is true) if the problem converges, and -1 otherwise.
"""
function solve_mpcc!(model, comp, plot_data = false, objective = :obj ;kwargs...) 
    # defining parameters
    gamma = get(kwargs, :gamma, 0.4);
    kappa = get(kwargs, :kappa, 0.2);
    sigma = get(kwargs, :sigma, 10);
    theta = get(kwargs, :theta, 10);
    mu = get(kwargs, :mu, 0.1);
    pi = get(kwargs, :pi, 1);
    tol = get(kwargs, :tol, 1e-8);
    k_max = get(kwargs, :k_max, 100);
    
    # getting model variables
    nc = length(comp);
    cvar = Array{Tuple{Any,Any},1}(undef, nc);
    cval = copy(cvar);
    nvc = Array{Int64,1}(undef, nc + 1);
    nvc[1] = 1;
    for i in 1:nc
        cvar[i] = (collect(Iterators.flatten(comp[i][1])), collect(Iterators.flatten(comp[i][2])));
        nvar = isa(cvar[i][1], VariableRef) ? 1 : length(cvar[i][1]); 
        nvc[i+1] = nvar + nvc[i];
    end
    obj = getindex(model, objective);

    # initial values
    k = 1;
    last_obj = 0;
    min_val = zeros(sum(nvc),1);
    ans = [];
    time = 0;
    if plot_data
        pl = [];
    end

    # setting objective function with penalty term for each pair of complementarity variables
    pen = Array{Any, 1}(undef, nc)
    for i in 1:nc
        if isa(cvar[i][1], VariableRef)
            pen[i] = @expression(model, cvar[i][1] * cvar[i][2]);
        else
            pen[i] = @expression(model, sum(cvar[i][1][j] * cvar[i][2][j] for j in 1:nvc[i+1]-nvc[i]));
        end
    end
    if isa(obj, NonlinearExpression)
        @NLobjective(model, Min, obj + pi * sum(pen[i] for i in 1:nc));
    else
        @objective(model, Min, obj + pi * sum(pen[i] for i in 1:nc));
    end

    # printing header for output
    @printf("\nIter   Termination_status    Objective    MPCC_compl lg(mu)   pi    CPUs \n")

    # outer loop
    while k <= k_max
        # reset variables for new iteration
        tol_c = mu^gamma;
        tol_p = theta*mu;
        set_optimizer_attributes(model, "mu_target" => mu, "mu_init" => mu,
                                        "dual_inf_tol" => tol_p, "constr_viol_tol" => tol_p,
                                        "compl_inf_tol" => tol_p, "print_level" => 0);
                                        # "compl_inf_tol" => tol_p);
        
        # solving the optimization problem                                
        optimize!(model);
        obj_val = objective_value(model);
        ans = value.(all_variables(model));
        time += solve_time(model);
        # println(value.(cvar[1][1]))
        for i in 1:nc
            cval[i] = (value.(cvar[i][1]), value.(cvar[i][2]));
            min_val[nvc[i]:nvc[i+1]-1] .= min.(cval[i][1], cval[i][2]);
        end
        compl = norm(min_val, Inf);

        # checking complementarity condition
        if compl <= tol_c
            @printf("%4.g %21s %.8e %.4e %2.2f %.1e, %.4f\n", k, termination_status(model), 
                    obj_val, compl, log10(mu), pi, time);
            if plot_data
                append!(pl, [[obj_val, mu, pi, compl]]);
            end
            # checking termination criteria
            if abs(obj_val - last_obj) < tol             
                println("\nIterations = ", k);
                println("Objective value = ", value(obj));     
                println("MPCC complementarity = ", compl);
                println("pi = ", pi);
                println("tol_p = ", tol_p);
                println("tol_c = ", tol_c);
                println("Time = ", time);
                break 
            else
                last_obj = obj_val;
            end
            k += 1;
            mu = mu*kappa;
            set_start_value.(all_variables(model), ans);
        else
            if pi < 1e14
                pi *= sigma;
                @objective(model, Min, obj + pi * sum(pen[i] for i in 1:nc));
            else
                error("Couldn't find a suitable value for pi")
            end
        end
    end
    if k > k_max
        println("Maximum number of outer iterations reached.");
        return -1
    end
    if plot_data
        return pl
    else
        return 0
    end
end







