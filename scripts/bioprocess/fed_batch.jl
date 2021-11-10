#
# Dynamic FBA for a fed batch process
# 
# Simple ecoli case study with 4 reactions and 3 external metabolites plus biomass
#

include("fed_batch_functions.jl");

### Input data
A = [0 -9.46 -9.84 -19.23;
     -35 -12.92 -12.73 0;
     -39.43 0 1.24 12.12];            # stoichiometry matrix
x0 = 0.001;
zglc_in = 5; 
nm, nr = size(A);
v0 = zeros(nr);              # TODO: estimate initial rates
tspan = [0; 11];
nfe = Int32(77);
dt = (tspan[2] - tspan[1])/nfe;

### Creating and solving the model
dfba = dfba_base_model(A, dt, ncp = 3);
objective_model_static(dfba);
# z0 = [10.8; 0.21; 0.];
z0 = [1.; 0.21; 0.];
D = dfba[:D];
fix.(D, 100.0/nfe/zglc_in, force = true);
# fix.(D, 0., force = true);
Z = Array{Float64,2}(undef, nfe+1, nm);
X = Array{Float64,1}(undef, nfe+1);
V = Array{Float64,2}(undef, nfe+1, nr);
Z[1,:] = z0;
X[1] = x0;
V[1,:] = v0;
for i in 1:nfe
    Z[i+1,:], X[i+1], V[i+1,:] = update_model_static(dfba, Z[i,:], X[i], V[i,:]);
    println("iteration $i: $(termination_status(dfba))");
end

# Plotting the results
fig = figure();
subplot(2,2,1);
for r in 1:nr
    plot(tspan[1]+dt:dt:tspan[2], V[2:end,r], label = "v$r");
end
legend();
subplot(2,2,2);
plot(tspan[1]:dt:tspan[2], Z[:,2], label = "O2");
legend();
subplot(2,2,3);
plot(tspan[1]:dt:tspan[2], Z[:,1], label = "Gluc");
plot(tspan[1]:dt:tspan[2], Z[:,3], label = "Acet");
legend();
subplot(2,2,4);
plot(tspan[1]:dt:tspan[2], X, label = "Biomass");
legend();

fig2 = figure()
plot(tspan[1]+dt:dt:tspan[2], sum(V, dims = 2)[2:end], label = "growth");
legend();