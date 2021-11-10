# MacMPEC Collection

This library contains some selected MPCC problems from the [MacMPEC Collection](https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC) originally built in AMPL by Sven Leyffer. For this library, these problems were written in Julia as a JuMP model and they are described below along with a link to the .jl file in the [github repository](https://github.com/carolinesnakama/mpcc_library/).

---
- ##### [bar-truss-3](https://github.com/Process-Optimization-and-Control/MPCCLibrary/tree/main/scripts/mpec_collection/bar-truss-3.jl)

Minimum weight design problem from M.C. Ferris and F. Tin-Loi, On the solution of a minimumWeight elastoplastic problem involving displacement and complementarity constraints, _Comp. Meth. in Appl. Mech & Engng_, 174:107-120, 1999.

Goal: minimize the volume of a struture with fixed topology that has to resist certain specified loads while keeping displacements within a specified limit. 
```math 
    \begin{aligned}
    \min \quad &L^Ta \\
    \text{s.t.} \quad &Q = SCw - SNz \\
                \quad &F = C^TQ \\
                \quad &w = -N^TQ + Hz + r \\
                \quad &S_i = Ea_i/L_i \text{ for } i = 1,\ldots,nm \\
                \quad &H_{i,jk} = 0.125S_j \text{ for } j = k \text{ and } i = 1,\ldots,nm \\
                \quad &H_{i,jk} = 0 \text{ for } j \neq k \text{ and } i = 1,\ldots,nm \\
                \quad &r_{i,j} = \sigma a_i \text{ for } j = 1,\ldots,ny \text{ and } i = 1,\ldots,nm \\
                \quad &Ta = 0 \\
                \quad &a \ge 0 \\
                \quad &-u^b \le u \le u^b \\
                \quad &0 \le w \perp z \ge 0
    \end{aligned}
```
where the variables are  

``u \in \mathbb{R}^{nd}`` - displacements  \
``z \in \mathbb{R}^{ny}`` - plastic multipliers  \
``a \in \mathbb{R}^{nm}`` - cross-sectional areas  \
``Q \in \mathbb{R}^{nm}`` - natural generalized stresses  \
``S \in \mathbb{R}^{nm}`` - element stiffness  \
``H_i \in \mathbb{R}^{ny \times ny}`` - hardening models  \
``r_i \in \mathbb{R}^{ny}`` - yield limits  \
``w_i(Q(z),z): \mathbb{R}^{ny} \rightarrow \mathbb{R}^{ny}`` - linear yield function  

and parameters are 

```math
    L = \begin{bmatrix} 500 & 400 & 500 \end{bmatrix}^T, \quad F = \begin{bmatrix} 400 & 600 \end{bmatrix}^T,
    \quad C = \begin{bmatrix} 0.6 & 0.0 & -0.6 \\ 0.8 & 1.0 & 0.8 \end{bmatrix}^T, \quad N = \begin{bmatrix} 1 & 1 & 1 \\ -1 & -1 & -1 \end{bmatrix}^T 
```
``E = 2e4``  \
``\sigma = 50``.  

``nm`` is the number of elements, ``nd`` is the number of structure degrees of freedom and ``ny`` is the number of yield functions per element. \
(``Ta = 0`` corresponds to technoogical constraints on the design variables ``a``)

---
- ##### [bard1](https://github.com/Process-Optimization-and-Control/MPCCLibrary/tree/main/scripts/mpec_collection/bard1.jl)

Two-level optimization toy model presented in J.F. Bard, Convex two-level optimization, _Mathematical Programming_ 40(1), 15-27, 1988.
```math
    \begin{aligned}
    \min_{x \ge 0} \quad &(x - 5)^2 + (2y + 1)^2 \\
    \text{s.t.} \quad &\min_{y \ge 0} \quad (y - 1)^2 - 1.5xy \\
    &\text{s.t.} \quad 3x - y \ge 3 \\
                     &\quad \quad -x + 0.5y \ge 4 \\
                     &\quad \quad -x -y \ge -7
    \end{aligned}
```
It can be reformulated as a MPCC problem by using the KKT conditions of the inner problem as constraints for the outer problem.
```math
    \begin{aligned}
    \min_{x,y \ge 0} \quad &(x - 5)^2 + (2y + 1)^2 \\
    \text{s.t.} \quad &2(y - 1) - 1.5x + \lambda_1 - 0.5\lambda_2 + \lambda_3 = 0 \\
                      &0 \le 3x - y -3 \perp \lambda_1 \ge  \\
                      &0 \le -x + 0.5y + 4 \perp \lambda_2 \ge 0 \\
                      &0 \le -x -y + 7 \perp \lambda_ 3 \ge 0
    \end{aligned}
```
where ``\lambda_i`` is the Lagrange multiplier corresponding to the ``i^\text{th}`` inequality of the inner problem.

---
- ##### [ralph2](https://github.com/Process-Optimization-and-Control/MPCCLibrary/tree/main/scripts/mpec_collection/ralph2.jl)

Toy model

```math
    \begin{aligned}
    \min_{x,y} \quad &x^2 + y^2 - 4xy \\
    \text{s.t.} \quad &0 \le x \perp y \ge 0
    \end{aligned}
```
---
- ##### [scale1](https://github.com/Process-Optimization-and-Control/MPCCLibrary/tree/main/scripts/mpec_collection/scale1.jl)

Toy model

```math
    \begin{aligned}
    \min_{x,y} \quad &(100x - 1)^2 + (y - 1)^2 \\
    \text{s.t.} \quad &0 \le x \perp y \ge 0
    \end{aligned}
```
---
