#
# Thermal energy storage problem
#

# using Revise;
using MPCCLibrary;
using PyPlot;

### Parameters

# Heat Parameters 
ρ_dh = 1000         #kg/m3
Cp_dh = 4.18        #kJ/(kg-K)
q_dh = 10           #m3/hr

# Temperature 
T_dh_ret = 30       #return
T_dh_minSup = 60    #minimum for supply
T_dh_maxSup = 70    #maximum for supply

# Capacity
V_tes = 100         #m3
V_whb = 10          #m3
V_phb = 10          #m3

# Bounds
Tlb = 30;
Tub = 100;

### Input
Q_whb = vcat(1.0*ones(10,1), 1.2*ones(10,1), 0.8*ones(10,1)) *1.254e6

### Defining sets for variables
S = ["T_tes", "T_phb", "T_whb"];         # state variables
C = ["q_whb", "q_A", "q_B", "Q_phb"];    # control variables
A = ["q_bp", "T_A", "T_B", "T_C"];       # aditional variables
 
### Discretization
ncp = 3;
tspan = (0.0, 30.0);
h = 1.0
nfe = Int32((tspan[2] - tspan[1])/h);
invM = collocation_mat(ncp);

### Defining model
model = create_model(linear_solver = "ma97");

### Variables
@variable(model, Tlb <= x[S, 1:nfe, 1:ncp] <= Tub, start = 60);
@variable(model, u[C, 1:nfe] >= 0, start = 0);
@variable(model, z[A, 1:nfe, 1:ncp] >= 0, start = 60);
@variable(model, x0[S]);

for i in A[2:end]
      set_lower_bound.(z[i,:,:], Tlb);
      set_upper_bound.(z[i,:,:], Tub);
end
set_lower_bound.(x["T_phb",:,:], T_dh_minSup);
set_upper_bound.(x["T_phb",:,:], T_dh_maxSup);
set_upper_bound.(u["q_whb",:], 15);
fix.(x0, 60, force=true);
set_start_value.(z["T_A",:,:],30.0);
set_start_value.(z["q_bp",:,:], 0.0);
set_start_value.(u["q_whb",:], 10.0);

### Mathematical model
tmp = Dict();
tmp["T_tes"] = @NLexpression(model, [i in 1:nfe, j in 1:ncp], u["q_A",i]/V_tes*(z["T_A",i,j] - x["T_tes",i,j])
                                                            + u["q_B",i]/V_tes*(z["T_B",i,j] - x["T_tes",i,j]));
tmp["T_phb"] = @NLexpression(model, [i in 1:nfe, j in 1:ncp], q_dh/V_phb*(z["T_C",i,j] - x["T_phb",i,j]) 
                                                             + u["Q_phb",i]/(V_phb*ρ_dh*Cp_dh));
tmp["T_whb"] = @NLexpression(model, [i in 1:nfe, j in 1:ncp], u["q_whb",i]/V_whb*(z["T_A",i,j] - x["T_whb",i,j]) 
                                                             + Q_whb[i]/(V_whb*ρ_dh*Cp_dh));
@NLexpression(model, ODE[k in S, i in 1:nfe, j in 1:ncp], tmp[k][i,j]);

### Model constraints
# discrtized ODE model
@NLconstraint(model, dODE0[k in S, j in 1:ncp], x[k,1,j] == x0[k] 
                                                          + h*(sum(invM[j,l]*ODE[k,1,l] for l in 1:ncp)));
@NLconstraint(model, dODE[k in S, i in 2:nfe, j in 1:ncp], x[k,i,j] == x[k,i-1,ncp] 
                                                                     + h*(sum(invM[j,l]*ODE[k,i,l] for l in 1:ncp)));
# algebraic constraints (balances)
@NLconstraint(model, mass[i in 1:nfe, j in 1:ncp], z["q_bp",i,j] == q_dh - u["q_A",i] + u["q_B",i] - u["q_whb",i]);
@NLconstraint(model, energ1[i in 1:nfe, j in 1:ncp], z["T_A",i,j] == (q_dh*T_dh_ret + u["q_B",i]*x["T_tes",i,j])/(q_dh + u["q_B",i]));
@NLconstraint(model, energ2[i in 1:nfe, j in 1:ncp], z["T_B",i,j] == (u["q_whb",i]*x["T_whb",i,j] + u["q_A",i]*x["T_tes",i,j])/(u["q_whb",i] + u["q_A",i]));
@NLconstraint(model, energ3[i in 1:nfe, j in 1:ncp], z["T_C",i,j] == (z["q_bp",i,j]*z["T_A",i,j] + (q_dh - z["q_bp",i,j])*z["T_B",i,j])/q_dh);

# Solving by simply adding a penalization term for the complementarity constraint
@NLobjective(model, Min, sum(1e-1*u["Q_phb",i]^2 for i in 1:nfe)                                     #goal: minimize fuel consumption
                         + sum(1e-2*u["q_A",i]*u["q_B",i] for i in 1:nfe)                            #complementarity penalty 
                        #  + 1e-4*sum((u[s,i] - u[s,i-1])^2 for s in ["q_A"; "q_B"] for i in 2:nfe)    #regularization for limiting big steps 
                         + 1e-4*sum((u[s,i] - u[s,i-1])^2 for s in ["q_A"] for i in 2:nfe)    #regularization for limiting big steps 
                         + 1e-5*sum((u[s,i] - u[s,i-1])^2 for s in ["q_B"] for i in 2:nfe)    #regularization for limiting big steps 
                         + 1e-3*sum((h*z["q_bp",i,j])^2 for i in 1:nfe for j in 1:ncp));             #regularization for preventing the usage of the bypass

optimize!(model);
cost = sum(value.(u["Q_phb",:])[i] for i in 1:nfe);

### Plotting results 
fig = figure(figsize = (15,15));

#Differential States
subplot(321);
ax1 = plot(tspan[1]:tspan[2], vcat(value.(x0["T_tes"]), value.(x["T_tes",:,end])), label = "T_tes");
ax1 = plot(tspan[1]:tspan[2], vcat(value.(x0["T_phb"]), value.(x["T_phb",:,end])), label = "T_phb");
ax1 = plot(tspan[1]:tspan[2], vcat(value.(x0["T_whb"]), value.(x["T_whb",:,end])), label = "T_whb");
title("State variables, @Obj $cost");
legend();

#Input data - heat load from waste heat boiler
subplot(322);
ax6 = step(tspan[1]:tspan[2]-1, Q_whb, label = "Q_whb (kJ/hr)")
title("Input data");
xlabel("t")
legend();

#Manipulated Variables
subplot(323);
ax2 = step(tspan[1]+1:tspan[2], value.(u["q_whb",:]), label = "q_whb")
ax2 = step(tspan[1]+1:tspan[2], value.(u["q_A",:]), label = "q_A")        
ax2 = step(tspan[1]+1:tspan[2], value.(u["q_B",:]), label = "q_B")
title("Control variables");
ylim([0.0, 15.0]);
legend();

subplot(324);
ax3 = step(tspan[1]+1:tspan[2], value.(u["Q_phb",:]), label ="Q_phb (kJ/hr)")
title("Heat from fuel");
ylim([0.0, 10.0]);
legend();

subplot(325);
#Algebraic States
ax4 = plot(tspan[1]+1:tspan[2], value.(z["q_bp",:,end]), label = "q_bypass")
title("Bypass flow");
ylim([0.0, 10.0]);
legend();

subplot(326);
ax5 = plot(tspan[1]+1:tspan[2], value.(z["T_A",:,end]), label = "T_A")  
ax5 = plot(tspan[1]+1:tspan[2], value.(z["T_B",:,end]), label = "T_B")  
ax5 = plot(tspan[1]+1:tspan[2], value.(z["T_C",:,end]), label = "T_C") 
title("Node temperatures")
legend();