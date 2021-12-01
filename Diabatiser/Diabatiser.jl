module Diabatiser
using QuadGK
using LinearAlgebra
using Optim

export lorentzian, diabatise, adiabatise, fit_diabat

function lorentzian(r, w, r0)
    return (1/(2*w))*(w^2/((r - r0)^2 + w^2))
end

<<<<<<< Updated upstream
function diabatise(r, a, f, p)
    phi, _ = quadgk(x -> f(x, p...), -Inf, r, rtol=1e-12, order=10)
    U = [
        cos(phi) -sin(phi);
        sin(phi)  cos(phi)
    ]
    d = adjoint(U)*a*U
    return d
=======
function nactomixingangle(r, f, p)
    phi, err = quadgk(x -> f(x, p...), -Inf, r, rtol=1e-15, order=10)
    return phi
>>>>>>> Stashed changes
end

@inline adiabatictodiabatic(phi) = [cos(phi) -sin(phi); sin(phi)  cos(phi)]

<<<<<<< Updated upstream
function adiabatise(r, d, f, p)
    phi, _ = quadgk(x -> f(x, p...), -Inf, r, rtol=1e-12, order=10)
    U = [
        cos(phi) -sin(phi);
        sin(phi)  cos(phi)
    ]
    a = U*d*adjoint(U)
    return a
end

=======
diabatise(a::Matrix, U::Matrix) = adjoint(U)*a*U
diabatise(a::Vector, U::Matrix) = adjoint(U)*a
adiabatise(a::Matrix, U::Matrix) = U*a*adjoint(U)
adiabatise(a::Vector, U::Matrix) = U*a

diabatise(a::T, phi::Float) where T<:Union{Matrix, Vector} = diabatise(a, adiabatictodiabatic(phi))
adiabatise(a::T, phi::Float) where T<:Union{Matrix, Vector} = adiabatise(a, adiabatictodiabatic(phi))

diabatise(a::T, r, f::Function, p) where T<:Union{Matrix, Vector} = diabatise(a, nactomixingangle(r, f, p))
adiabatise(a::T, r, f::Function, p) where T<:Union{Matrix, Vector} = adiabatise(a, nactomixingangle(r, f, p))

>>>>>>> Stashed changes
@inline so_derivative(r, a) = 2*(diff(diff(a) ./ diff(r))) ./ (diff(r)[1:length(r)-2] .+ diff(r)[2:length(r)-1])

function get_loss(r, a, f, p)
    d = map(i -> diabatise(a[i, :, :], r[i], f, p), collect(1:length(r)))
    loss = sum(abs.(so_derivative(r, [d_[1, 1] for d_ in d])))
    return loss
end

function detect_r0(r, a)
    so_grads = so_derivative(r, a)
    ind_max = findmax(abs.(so_grads))[2]
    return r[ind_max + 1]
end

function fit_diabat(a::Array, r, f; w0=0.01)
    r0 = map(i -> detect_r0(r, a[:, i, i]), [1, 2])
    r0 = sum(r0)/length(r0)
    o = optimize(p -> get_loss(r, a, f, p), [w0, r0])
    return Optim.minimizer(o)
end

function fit_diabat(r, a::Vector, f; w0=0.01)
    adiabats = Array{Float64}(undef, length(a), 2, 2)
    for i=1:length(r)
        adiabats[i, :, :] = a[i]
    end
    return fit_diabat(r, adiabats, f; w0=w0)
end

end