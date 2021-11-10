#  
#  Flash tank problem
#

# using Revise;
using MPCCLibrary;
using PyPlot;

### Parameters

# Antoine's equation
A = [3.97786; 4.00139; 3.93002];
B = [1064.84; 1170.874; 1182.774];
C = [-41.136; -48.833; -52.532];

### Input data

F = 1;                   #feed
z = [0.5; 0.3; 0.2];     #feed concentration
Tt = 380:400;            #temperature range
p = 5;                   #bar
nc = length(z);
nt = length(Tt);

### Defining model
model = create_model(linear_solver = "ma27", print_level = 5);

# model variables
@variable(model, 0 <= x[1:nc, 1:nt] <= 1);
@variable(model, 0 <= y[1:nc, 1:nt] <= 1);
@variable(model, L[1:nt] >= 0, start = .5);
@variable(model, V[1:nt] >= 0, start = .5);
@variable(model, T[1:nt] >= 0);
@variable(model, K[1:nc, 1:nt] >= 0);
@variable(model, psat[1:nc, 1:nt] >= 0);      #pure component vapor pressure
# auxiliary variables
@variable(model, log_psat[1:nc, 1:nt] >= 0);
@variable(model, 0 <= a[1:nt] <= 1);          #corrected vapor fraction
@variable(model, at[1:nt]);                   #vapor fration calculated with the Rachford-rice equation
@variable(model, sv[1:nt] >= 0);              #slack for vapor fraction
@variable(model, sl[1:nt] >= 0);              #slack for liquid fraction

#initial guesses
set_start_value.(x, z);
set_start_value.(y, z);
set_start_value.(T, Tt);

## constraints
@NLconstraint(model, ant[i in 1:nc, j in 1:nt], log_psat[i,j] == A[i] - B[i]/(T[j] + C[i]));                #Antoine's equation
@NLconstraint(model, pv[i in 1:nc, j in 1:nt], psat[i,j] == 10^log_psat[i,j]);
@NLconstraint(model, eqp[i in 1:nc, j in 1:nt], K[i,j] == psat[i,j]/p);                                     #Raoult's law
@NLconstraint(model, eqc[i in 1:nc, j in 1:nt], K[i,j]*x[i,j] == y[i,j]);
@NLconstraint(model, balc[i in 1:nc, j in 1:nt], L[j]*x[i,j] + V[j]*y[i,j] == F*z[i]);                      #Balance per component
@NLconstraint(model, balt[j in 1:nt], L[j] + V[j] == F);                                                    #Total balance
@NLconstraint(model, rr[j in 1:nt], sum((z[i]*(K[i,j] - 1))/(1 + at[j]*(K[i,j] - 1)) for i in 1:nc) == 0);  #Rachford-Rice equation
@constraint(model, aux[j in 1:nt], a[j] - sv[j] + sl[j] - at[j] == 0);                                      
@constraint(model, temp[j in 1:nt], T[j] == Tt[j]);                                                         #Temperature profile

## objective function defined as an expression
comp = [(sv, V); (sl, L)];                                                 #complementarity pairs
@expression(model, obj, sum(0.5*(a[j]*F - V[j])^2 for j in 1:nt));
solve_mpcc!(model, comp);

###Plotting the results
fig = figure();

subplot(211);
ax1 = plot(Tt, value.(a), label = "V/F corrected", marker = "o");
ax1 = plot(Tt, value.(at), label = "V/F", marker = "o");
xlabel("T (K)");
legend();

subplot(212);
ax2 = plot(Tt, value.(sv), label = "sV", marker = "o");
ax2 = plot(Tt, value.(sl), label = "sL", marker = "o");
xlabel("T (K)");
legend();




