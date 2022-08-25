# --- W. Somogyi ---

module Diabatiser
using QuadGK
using LinearAlgebra
using Optim

export lorentzian, diabatise, adiabatise, fit_diabat

"""
Example Lorentzian NAC function parameterised by width `w` and central geometry `rc`.
    
Obeys the normalisation condition  ∫f(r) dr = π/2
"""
function lorentzian(r, w, rc)
    return (1/2)*((w/2)/((r - rc)^2 + (w/2)^2))
end

"""
Integrate a NAC to obtain the mixing angle at a given geometry

# Arguments
- `r`: the geometry at which to compute the mixing angle
- `f`: the function describing the NAC
- `p`: vector of parameters for the function `f`
"""
function nactomixingangle(r, f, p)
    phi, err = quadgk(x -> f(x, p...), -Inf, r, rtol=1e-12, order=10)
    return phi
end

"""
Create the adiabatic to diabatic transformation matrix from a mixing angle `phi`
"""
@inline adiabatictodiabatic(phi) = [cos(phi) -sin(phi); sin(phi) cos(phi)]

"""
Diabatise a matrix or vector `a` according to the transformation matrix `U`
"""
diabatise(a::Matrix, U::Matrix) = adjoint(U)*a*U
diabatise(a::Vector, U::Matrix) = adjoint(U)*a

"""
Diabatise a matrix or vector `a` according to the NAC term `phi`.
"""
diabatise(a::T, phi::Real) where T<:Union{Matrix, Vector} = diabatise(a, adiabatictodiabatic(phi))

"""
Diabatise a matrix or vector according to a NAC function.

# Arguments
- `a`::Union{Matrix, Vector}: the adiabatic vector or matrix at geometry `r`
- `r`: the geometry at which to perform the didabatisation
- `f`::callable: the function describing the NAC as a function of geometry `r`
- `p`: the parameters of the function `f`
"""
diabatise(a::T, r, f::Function, p) where T<:Union{Matrix, Vector} = diabatise(a, nactomixingangle(r, f, p))

"""
Adiabatise a matrix or vector - the inverse of `diabatise`.

See [`diabatise!`](@ref)
"""
adiabatise(a::Matrix, U::Matrix) = U*a*adjoint(U)
adiabatise(a::Vector, U::Matrix) = U*a
adiabatise(a::T, phi::Real) where T<:Union{Matrix, Vector} = adiabatise(a, adiabatictodiabatic(phi))
adiabatise(a::T, r, f::Function, p) where T<:Union{Matrix, Vector} = adiabatise(a, nactomixingangle(r, f, p))

"""
Calculate the second derivative of a vector `a` over a grid of geometries `r`.
"""
@inline so_derivative(a, r) = 2*(diff(diff(a) ./ diff(r))) ./ (diff(r)[1:length(r)-2] .+ diff(r)[2:length(r)-1])

"""
Calculate the sum of second derivatives of the diabatic terms given a set of adiabatic terms.

# Arguments
- `a`: the adiabatic terms on a grid of geometries `r`
- `r`: the grid of geometries
- `f`: the function describing the NAC
- `p`: the variables parameterising the NAC to pass as arguments to `f`
"""
function get_loss(a, r, f, p)
    d = map(i -> diabatise(a[i, :, :], r[i], f, p), collect(1:length(r)))
    loss = sum(abs.(so_derivative([d_[1, 1] for d_ in d], r)))
    return loss
end

"""
Convenience function that attempts to make an initial guess for the central geometry.
"""
function detect_rc(a, r)
    so_grads = so_derivative(a, r)
    ind_max = findmax(abs.(so_grads))[2]
    return r[ind_max + 1]
end

"""
Top-level function to fit a NAC function and return the fitted parameters

# Arguments
-`a`: the adiabatic terms on a grid of geometries `r`
-`r`: the grid of geometries
-`f`: the function describing the NAC
-`w0=0.01`: initial guess for the characteristic width of the NAC
-`rc=nothing`: initial guess for the central geometry, if `nothing` then an attempt is made to determine it automatically
"""
function fit_diabat(a::Array, r, f; w0=0.01, rc0=nothing)
    rc0 === nothing ? rc0 = detect_rc(a[:, 1, 1], r) : rc0
    o = optimize(p -> get_loss(a, r, f, p), [w0, rc0])
    return Optim.minimizer(o)
end

function fit_diabat(a::Vector, r, f; w0=0.01, rc0=nothing)
    adiabats = Array{Float64}(undef, length(a), 2, 2)
    for i=1:length(r)
        adiabats[i, :, :] = a[i]
    end
    return fit_diabat(adiabats, r, f; w0=w0, rc0=rc0)
end

end
