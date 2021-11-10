#
# MacMPEC Collection by Sven Leyffer 
# Problem: ralph2 
#

# using Revise;
using MPCCLibrary;
using LaTeXStrings, PyPlot;   # used for plotting the algorithm's progress

model = create_model(print_level = 0);
@variable(model, x >= 0, start = 1);
@variable(model, y >= 0, start = 1);

# define original objective function as an expression
@expression(model, obj, x^2 + y^2 - 4 * x * y);     

ans = solve_mpcc!(model, [(x, y)], true);

# plotting the results
n = length(ans)
nd = length(ans[1])
data = Array{Float64,2}(undef, n, nd)
for i in 1:n
    data[i,:] = ans[i]
end
fig = figure();
ax = PyPlot.subplot();
PyPlot.xlabel("Outer iterations");
PyPlot.ylabel(L"log_{10}");
for i in 1:nd
    PyPlot.plot(1:n, log.(data[:,i]), marker = i+6);
end
PyPlot.grid("on")
PyPlot.legend(["Objective", L"μ_k", L"π_k", "max(complementarity)"]);

