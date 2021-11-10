#
# Bioprocess Optimization Problem (dFBA as inner problem)
# 
# Simple ecoli case study with 4 reactions and 3 external metabolites plus biomass
#

include("fed_batch_functions.jl");

### Input data
A = [0 -9.46 -9.84 -19.23;
     -35 -12.92 -12.73 0;
     -39.43 0 1.24 12.12];            # stoichiometry matrix
x0 = 0.001;
# z0 = [10.8; 0.21; 0.4];
z0 = [1.; 0.21; 0.4];
zglc_in = 5; 
nm, nr = size(A);
v0 = zeros(nr);              # TODO: estimate initial rates
tspan = [0; 11];
nfe = Int32(77);
dt = (tspan[2] - tspan[1])/nfe;

### Creating the inner model
dfba = dfba_base_model(A, dt, nfe = nfe, print_level = 5);
x = dfba[:x];
z = dfba[:z];
v = dfba[:v];
D = dfba[:D];
all_vars = all_variables(dfba)
dvars = all_vars[1:length(all_vars)-length(D)];  # min w.r.t. x, z, v (not D)
fix(x[1,1], x0, force = true);
fix.(v[1,:], v0, force = true);
fix.(z[1,:,1], z0, force = true);
objective_model_static(dfba);

### Setting up the outter model
comp = set_bilevel_opt_problem!(dfba, dvars);
# @objective(dfba, Min, -x[end,end]);
@expression(dfba, obj, -x[end,end]);
@constraint(dfba, lim_D, sum(D[i] for i in 1:nfe) == 100.0/zglc_in);
@constraint(dfba, D_lb[i in 1:nfe], D[i] >= 0);
#solving with Leyffer's algorithm
solve_mpcc!(dfba, comp, tol = 1e-6);

#solving by adding constraints for the complementarity 
# for c in comp
#      @constraint(dfba, c[1]*c[2] <= 1e-5);
# end
# optimize!(dfba);

# Plotting the results
fig = figure();
subplot(2,2,1);
for r in 1:nr
    plot(tspan[1]+2*dt:dt:tspan[2], value.(v[2:end,r]), label = "v$r");
end
legend();
subplot(2,2,2);
plot(tspan[1]:dt:tspan[2], vcat(z0[2], value.(z[:,2,end])), label = "O2");
legend();
subplot(2,2,3);
plot(tspan[1]:dt:tspan[2], vcat(z0[1], value.(z[:,1,end])), label = "Gluc");
plot(tspan[1]:dt:tspan[2], vcat(z0[3], value.(z[:,3,end])), label = "Acet");
legend();
subplot(2,2,4);
plot(tspan[1]:dt:tspan[2], vcat(x0, value.(x[:,end])), label = "Biomass");
legend();

fig2 = figure()
plot(tspan[1]+2*dt:dt:tspan[2], sum(value.(v), dims = 2)[2:end], label = "growth");
legend();