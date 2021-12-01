module Diabatiser
using QuadGK
using LinearAlgebra
using Optim

export lorentzian, diabatise, adiabatise, fit_diabat

function lorentzian(r, w, r0)
    return (1/(pi*w))*(w^2/((r - r0)^2 + w^2))
end

function nactomixingangle(r, f, p)
    phi, err = quadgk(x -> f(x, p...), -Inf, r, rtol=1e-12, order=10)
    return phi
end

@inline adiabatictodiabatic(phi) = [cos(phi*pi/2) -sin(phi*pi/2); sin(phi*pi/2) cos(phi*pi/2)]

diabatise(a::Matrix, U::Matrix) = adjoint(U)*a*U
diabatise(a::Vector, U::Matrix) = adjoint(U)*a
adiabatise(a::Matrix, U::Matrix) = U*a*adjoint(U)
adiabatise(a::Vector, U::Matrix) = U*a

diabatise(a::T, phi::Real) where T<:Union{Matrix, Vector} = diabatise(a, adiabatictodiabatic(phi))
adiabatise(a::T, phi::Real) where T<:Union{Matrix, Vector} = adiabatise(a, adiabatictodiabatic(phi))

diabatise(a::T, r, f::Function, p) where T<:Union{Matrix, Vector} = diabatise(a, nactomixingangle(r, f, p))
adiabatise(a::T, r, f::Function, p) where T<:Union{Matrix, Vector} = adiabatise(a, nactomixingangle(r, f, p))

@inline so_derivative(a, r) = 2*(diff(diff(a) ./ diff(r))) ./ (diff(r)[1:length(r)-2] .+ diff(r)[2:length(r)-1])

function get_loss(a, r, f, p)
    d = map(i -> diabatise(a[i, :, :], r[i], f, p), collect(1:length(r)))
    loss = sum(abs.(so_derivative(r, [d_[1, 1] for d_ in d])))
    return loss
end

function detect_r0(a, r)
    so_grads = so_derivative(a, r)
    ind_max = findmax(abs.(so_grads))[2]
    return r[ind_max + 1]
end

function fit_diabat(a::Array, r, f; w0=0.01, r0=nothing)
    r0 == nothing ? r0 = detect_r0(a[:, 1, 1], r) : r0
    o = optimize(p -> get_loss(a, r, f, p), [w0, r0])
    return Optim.minimizer(o)
end

function fit_diabat(a::Vector, r, f; w0=0.01, r0=nothing)
    adiabats = Array{Float64}(undef, length(a), 2, 2)
    for i=1:length(r)
        adiabats[i, :, :] = a[i]
    end
    return fit_diabat(adiabats, r, f; w0=w0, r0=r0)
end

end