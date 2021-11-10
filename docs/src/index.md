# Introduction

Mathematical programs with complementarity constraints (MPCC) are optimization problems of the form
```math
    \begin{aligned}
    \min \quad &\varphi := f(x, y, z) \\
    \text{s.t.} \quad &h(x, y, z) = 0 \\
    &g(x, y, z) \ge 0 \\
    &0 \le x \perp y \ge 0
    \end{aligned}  
```
in which the last constraint means that either ``x \in \mathbb{R}^n \ge 0`` or ``y \in \mathbb{R}^n \ge 0``. They can be used to model processes that, for example, present switches or nonsmooth decisions with the advantage of incorporating discrete decisions in a single nonlinear programming (NLP) formulation[^1]. 

[^1]:  Lorenz T. Biegler. _Nonlinear programming: concepts, algorithms, and applications to chemical processes_. Vol. 10. Siam, 2010

This library is intended to be a compilation of MPCC (Mathematical Programs with Complementarity Constraints) problems to test algorithms developed in [Julia](https://julialang.org) using [JuMP](https://jump.dev/JuMP.jl/stable/) as modeling language and to be used as case studies. It also contains an implementation of the first algorithm proposed by Leyffer et al (2006)[^2] that can be used to solve examples in this library.

[^2]: [Leyffer, S., López-Calva, G. and Nocedal, J. Interior methods for mathematical programs with complementarity constraints. _SIAM Journal on Optimization_, 17(1), pp.52-77, 2006.](https://epubs.siam.org/doi/abs/10.1137/040621065)

```@contents 
Pages = [
    "pages/flash_problem.md"
    "pages/tes_problem.md"
    "pages/bio_problem.md"
    "pages/biopt_problem.md"
    "pages/mpec_collection.md"
    "pages/mpcc_methods.md"
    ]
Depth = 4
```
