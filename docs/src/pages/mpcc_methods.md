# Methods

### Dealing with MPCCs

```@docs
create_model
solve_mpcc!
```

### Auxiliary methods

#### Orthogonal collocation

Discretization method for handling dynamic optimization.

```@docs
collocation_mat
collocation_matrix
```

#### Bilevel optimization

To formulate a bilevel optimization problem as a MPCC problem, we need to add the first order conditions of the inner problem as constraints of the outer optimization problem. 

```@docs
set_bilevel_opt_problem!
```