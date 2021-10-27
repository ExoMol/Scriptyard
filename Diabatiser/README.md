# Diabatiser
_by Wilfrid Somogyi_

This is a Julia module for transforming potential energy curves, dipole moment curves, spin-orbit coupling curves, etc. between diabatic and adiabatic representations. It provides methods for fitting arbitrary functional forms of the non-adiabatic coupling curves in order to obtain the unitary transformation matrix between the two representations.

A concise discussion of the adiabatic-diabatic transformation can be found in [_Abrol and Kuppermann_, J. Chem. Phys. 116, 1035-1062 (2002)](https://doi.org/10.1063/1.1419257). Their discussion relates to solving the Poisson equation to obtain a multi-dimensional NAC for the triatomic system H<sub>3</sub>, but the material in section I. is independent of the molecule and applies generally to any two-state system.

## Setup

Download the appropriate Julia language executable from https://julialang.org/downloads/ and add the executable to your PATH variable if neccessary (usually only required for Linux)

```
export PATH="$PATH:/path/to/julia/bin"
```

Run Julia (type `julia` in the command line) and install the required packages ([`Optim.jl`](https://github.com/JuliaNLSolvers/Optim.jl) and [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl)) by pressing `]` in the Julia REPL to access the built-in package manager and adding the packages

```
julia> ]
(@v1.6) pkg> add Optim QuadGK
```

You may also want to install the packages `CSV.jl`, `DataFrames.jl` and `LinearAlgebra.jl`. While not scrictly necessary, they are useful for reading and pre-processing data from CSV files. The `matplotlib.pyplot` binding module `PyPlot.jl` is also required to run the `test.jl` file included in the repository.

To make use of the `Diabatiser` module in your own Julia script, add the following import statement to the top of your `.jl` file with the path to the location of `Diabatiser.jl`, e.g

```julia
include("/path/to/Scriptyard/Diabatiser/Diabatiser.jl")
```

## Usage

### Diabatising Adiabatic Curves

To diabatise the potential energy operator for a given geometry we can call the function `Diabatiser.diabatise`, for example

```julia
diabat = Diabatiser.diabatise(r, adiabat, nac_func, p)
```

Which performs the diabatisation at internuclear distance `r` of the 2x2 potential energy operator `adiabat` parametrised by the non-adiabatic coupling function `nac_func` with parameters specified by the list `p`. A more concrete example is given below showing the diabatisation of the potential energy operator at a geometry (r = 1.580) close to the avoided crossing paramterised by a Lorentzian (width = 0.012Å, r0 = 1.588Å) NAC function.

```julia
adiabat = [-149.791 0; 0 -149.784] #2x2 matrix with adiabatic potential energies on the diagonal
p = [0.012, 1.588] #width and central value parameters for Lorentzian function
diabat = Diabatiser.diabatise(1.580, adiabat, Diabatiser.lorentzian, p)
```

Which returns a 2x2 matrix with the diabatic potential energies and the off-diagonal diabatic couplings (V1 = -149.789, V2 = -149.786, V12 = V21 = 0.00291218). If the adiabatic potential energy matrices for each geometry are stored as a vector of matrices, the function can be vectorised using standard Julia dot notation

```julia
adiabats = [[-149.791 0; 0 -149.782], [-149.791 0; 0 -149.784], [-149.792 0; 0 -149.786]]
rs = [1.575, 1.580, 1.589]
p = [0.012, 1.588]
diabats = Diabatiser.diabatise.(rs, adiabats, Diabatiser.lorentzian, [p for i=1:length(adiabats)])
```

### Adiabatising Diabatic Curves

In principle adiabatisation can be performed with no additional knowledge of the non-adiabatic couplings by simply diagonalising the diabatic potential energy matrix. However, `Diabatiser.jl` provides a function for explicit adiabatisation using the inverse transformation obtained from the NAC function. The syntax is identical to the previous example, except the adiabatic and diabatic matrices are swapped

```julia
adiabat = Diabatiser.adiabatise(r, diabat, nac_func, p)
```

### Optimising the NAC Function

In the case that the parameters of the NAC are not known prior to the diabatisation, `Diabatiser.jl` provides the function `fit_diabat` that minimises a loss function to obtain the best-fit parameters for an arbitrary functional form of the NAC function. In this case, the entire vector of geometries is processed simultaneously, rather than on a geometry-by-geometry basis.

```julia
opt = Diabatiser.fit_diabat(rs, adiabats, Diabatiser.lorentzian)
```

The adiabats can alternatively be provided as a 3-dimensional array, with geometries along the first dimension. 

The loss function is the sum of second order derivatives of the PEC, which aims to minimise gradient changes in the region of the avoided crossing. Because the parameter space of the NAC often results in numerous local minima in the loss function, particularly when the position of the avoided crossing is not well constrained, the optimisation procedure automatically guesses an initial centre for the avoided crossing by searching for the geometry corresponding to the largest second order derivative. Providing an initial guess for the width is more complicated, but when the position of the avoided crossing is well constrained by an initial guess, the loss function generally has only a few local minima and an initial guess for the width of the NAC with the correct order of magnitude is usually sufficient. Thus the initial guess for the width of the NAC defaults to `w=0.01`, which generally gives correct minimisation for bond lengths in Angstrom or Bohrs, an alternative initial guess can be provided via the keyword argument `w`, this is useful if the length dimensions are much larger or smaller than Angstroms e.g nanometres

```julia
opt = Diabatiser.fit_diabat(rs, adiabats, Diabatiser.lorentzian; w0=0.001)
```

### Functional Forms of the NAC

Although the Lorentzian function (`Diabatiser.lorentzian`) supplied as part of the module is usually a good approximation of the NAC, any arbitrary user-defined functional form parametrised by a HWHM (`w`) and central value (`r0`) can in principle be provided to the (a)diabatisation and fit functions. This could be, for example, a Gaussian distribution.


