# Julia MPCC Library

Library of MPCC (Mathematical Programs with Complementarity Constraints) problems written in [Julia](https://julialang.org) using [JuMP](https://jump.dev/JuMP.jl/stable/). 

This library also contains an implementation of the first algorithm in [Leyffer et al (2006)](https://epubs.siam.org/doi/abs/10.1137/040621065) to solve MPCCs with an interior point solver. 

## Documentation

Documentation includes a description of the problems and a reference to the available functions. You can find it [here](https://carolinesnakama.github.io)
or build it yourself. For the latter, you can download this repository and generate the documentation with the source files available. 

1- Make sure you have the Documenter package. Press `]` in the julia REPL to start the Pkg REPL
```julia
(@v1.5) pkg> add Documenter
```

2- Run the `make.jl` file in the `docs` folder
```julia
julia> include("/path/to/MPCCLibrary/docs/make.jl")
```

3- The pages will be available on `/path/to/MPCCLibrary/docs/build/`. Just open `index.html` to navigate through the documentation.

## How to use this library

You can just download the problem you want from the `scripts` folder. The flash tank problem and the ones from the MacMPEC Collection were written using the algorithm from Leyffer et al (2006) to solve the problem. You can modify them and customize your solver and objective function or run the file as it is. For the latter, you need to add this library as a package in your Julia installation. To do so, you can follow these steps:

1- Download this repository.

2- Start julia and press `]` to enter the Pkg REPL.

3- Add this package 
```julia
(@v1.5) pkg> add /path/to/MPCCLibrary/
```

4- Press `backspace` to return to the julia REPL.

5- Run the script you wish.

## Problems
- Flash Tank Problem
- Thermal Energy Storage Problem
- Bioprocess Optimization Problem
- Bilevel Optimization Problem
- [MacMPEC Collection](https://wiki.mcs.anl.gov/leyffer/index.php/MacMPEC) by Sven Leyffer (originally written in AMPL).
