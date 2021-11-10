# Flash Tank Problem

#### Description
For this example[^1], a flash tank wih feed flow ``F`` and composition ``z_i`` with ``i = \mathcal{C}`` is considered. We want to analyze how the split between vapor ``V`` and liquid ``L`` products with composition ``y_i`` and ``x_i`` respectively vary with temperature for a fixed pressure ``p``. We can use the Rachford-Rice equation 

[^1]: [Kungurtsev, V. and Jäschke, J. Pathfollowing for Parametric Mathematical Programs with Complementarity Constraints. _ePrints for the Optimization Community_, 2019.](http://www.optimization-online.org/DB_HTML/2019/02/7084.html)

```math
    \sum_{i \in \mathcal{C}}\frac{z_i(K_i - 1)}{1+a_t(K_i - 1)} = 0
```
to calculate the fraction of the feed that goes to the vapor phase (represented by ``a_t``) with ``K_i`` given by Raoult's law
```math
    K_i = \frac{p_i^\text{sat}(T)}{p} = \frac{y_i}{x_i} \quad \text{for }i \in \mathcal{C}.
```
``p_i^\text{sat}(T)`` is the vapor pressure of the pure component ``i`` at temperature ``T`` calculated using Antoine's equation
```math
    \log_{10}(p_i^\text{sat}) = A_i - \frac{B_i}{T + C_i}
```
where ``A_i``, ``B_i`` and ``C_i`` are constants for each compound ``i``. 

A problem with this formulation is that, if T is lower than the mixture's bubble point or larger than its dew point, ``a_t`` assumes negative values and values greater than one respectively and, since ``a_t`` represents the ratio ``V/F``, that would be physically impossible. We can deal with it by formulating this problem as an optimization problem with complementarity constraints that enforce ``V/F`` to be ``0`` if ``T`` is lower than bubble point and ``1`` if it is larger than dew point. 

The optimization problem is then modeled as 
```math
    \begin{aligned}
        \min \quad &\frac{1}{2}(aF - V)^2 dt \\
                s.t. \quad &\sum_{i \in \mathcal{C}}\frac{z_i(K_i - 1)}{1+a_t(K_i - 1)} = 0 \\
                           &K_i = \frac{p_i^\text{sat}}{p} \quad \text{for }i \in \mathcal{C} \\
                           &K_i = \frac{y_i}{x_i} \quad \text{for }i \in \mathcal{C} \\
                           &\log_{10}(p_i^\text{sat}) = A_i - \frac{B_i}{T + C_i} \quad \text{for }i \in \mathcal{C} \\
                           &L + V = F \\
                           &Lx_i + Vy_i = Fz_i \quad \text{for }i \in \mathcal{C} \\
                           &a - s_V + s_L - a_t = 0 \\
                           &0 \le s_V \perp V \ge 0 \\
                           &0 \le s_L \perp L \ge 0 \\
                           &0 \le a \le 1 \\
                           &L, V \ge 0 \\
                           &K_i, p_i^\text{sat}, x_i, y_i \ge 0 \quad \text{for }i \in \mathcal{C}.
    \end{aligned}
```
The 5``^\text{th}`` and 6``^\text{th}`` constraints correspond to material balances, total and componentwise respectively. Constraints 7-10 describe the complementarity conditions, ``s_L`` and ``s_V`` are slacks variables that represent how much ``a_t`` is lower than ``0`` and larger than ``1`` respectively. Variable ``a``, then, represents the actual ratio ``V/F``, which is enforced by the 7``^\text{th}`` and 10 ``^\text{th}`` constraints. Note that ``a = V/F`` is not enforced as a hard constraint and, instead, is used as the objective function. 

#### Implementation
The implementation in [flash.jl](https://github.com/Process-Optimization-and-Control/MPCCLibrary/blob/main/scripts/flash.jl) uses

``F = 1``, \
``T \in [380,400]`` in increments of ``1`` K, \
and 3 compounds with the following Antoine's constants and inlet composition

| Compound | A       | B        | C       | z   |
|:--------:|:-------:|:--------:|:-------:|:---:|
| 1        | 3.97786 | 1064.84  | -41.136 | 0.5 |
| 2        | 4.00139 | 1170.875 | -48.833 | 0.3 |
| 3        | 3.93002 | 1182.774 | -52.532 | 0.2 |

This problem is solved using the `solve_mpcc!` function provided with the default parameters.  




