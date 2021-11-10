# Bilevel Optimization Problem

Two toy models[^1] are implemented as simple examples of bilevel optimization problems. They are solved by being reformulated as MPCCs problems, in which the first-order conditions of the inner problem are added as constraints of the outer problem. Implementation of the following models can be found in [bilevel.jl](https://github.com/Process-Optimization-and-Control/MPCCLibrary/tree/main/scripts/bilevel.jl)

[^1]: [Dempe, S. Bilevel optimization: Reformulation and first optimality conditions. _Generalized Nash Equilibrium Problems, Bilevel Programming and MPEC_, pp. 1-20. Springer, 2017.](https://link.springer.com/chapter/10.1007/978-981-10-4774-9_1)

##### Toy Model 1

```math
    \begin{aligned}
    \min_{x,y} \quad &(x - 1)^2 + (y - 1)^2 \\
    \text{s.t.} \quad &\min_{y} \quad -y \\
    &\text{ s.t.} \quad x + y \le 1 \\
    &\quad \quad - x + y \le 1
    \end{aligned}
```

##### Toy Model 2

```math
    \begin{aligned}
    \min_{a,b} \quad &- a - 2b \\
    \text{s.t.} \quad &2a - 3b \ge -12 \\
    &a + b \le 14 \\
    &\min_{b} \quad -b \\
    &\text{ s.t.} \quad -3a + b \le -3 \\
    &\quad \quad \quad 3a + b \le 30
    \end{aligned}
```