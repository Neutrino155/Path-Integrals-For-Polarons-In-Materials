# Feynman.jl

function polaron_propagator(τ, v, w)
    c = (v^2 - w^2) / v^3
    result = (1 - v * c) * τ + c * (phonon_propagator(0, v) - phonon_propagator(τ, v))
    return isa(result, Complex) ? result : result > 0 ? result : eps()
end

function polaron_propagator(τ, v, w, β)
    c = (v^2 - w^2) / v^3
    result = c * (phonon_propagator(0, v, β) - phonon_propagator(τ, v, β)) + (1 - c * v) * τ * (1 - τ / β)
    return isa(result, Complex) ? result : result > 0 ? result : eps()
end

function polaron_propagator(τ, v::Vector, w::Vector, β)
    result =  τ * (1 - τ / β) + sum((h_i(i, v, w) / v[i]^2) * ((1 + exp(-v[i] * β) - exp(-v[i] * τ) - exp(v[i] * (τ - β))) / (v[i] * (1 - exp(-v[i] * β))) - τ * (1 - τ / β)) for i in eachindex(v))
    return isa(result, Complex) ? result : result > 0 ? result : eps()
end

function polaron_propagator(τ, v::Vector, w::Vector)
    result =  τ + sum((h_i(i, v, w) / v[i]^2) * ((1 - exp(-v[i] * τ)) / v[i] - τ) for i in eachindex(v))
    return isa(result, Complex) ? result : result > 0 ? result : eps()
end

S₀(v, w) = -(3 / (4 * v)) * (v^2 - w^2)

S₀(v::Vector, w::Vector) = -3 * sum(C_ij(i, j, v, w) / (v[j] * w[i]) for i in eachindex(v), j in eachindex(w))

S₀(v, w, β) = -3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β)) 

S₀(v::Vector, w::Vector, β) = -3 * sum(C_ij(i, j, v, w) / (v[j] * w[i]) * (coth(β * v[j] / 2) - 2 / (β * v[j])) for i in eachindex(v), j in eachindex(w))

E₀(v, w) = 3 * (v - w) / 2

E₀(v::Vector, w::Vector) = 3 * sum(v .- w) / 2

E₀(v, w, β) = 3 * (v - w) / 2 - 3 / β * (log(v / w) - log(1 - exp(-v * β)) + log(1 - exp(-w * β))) 

E₀(v::Vector, w::Vector, β) = 3 * sum(v .- w) / 2 - 3 / β * sum(log(v[i] / w[i]) - log(1 - exp(-v[i] * β)) + log(1 - exp(-w[i] * β)) for i in eachindex(v))

function κ_i(i, v::Vector, w::Vector)
    κ = v[i]^2 - w[i]^2
    κ *= prod(j != i ? (v[j]^2 - w[i]^2) / (w[j]^2 - w[i]^2) : 1 for j in eachindex(v))
    return κ
end

function h_i(i, v::Vector, w::Vector)
    h = v[i]^2 - w[i]^2
    h *= prod(j != i ? (w[j]^2 - v[i]^2) / (v[j]^2 - v[i]^2) : 1 for j in eachindex(v))
    return h
end

function C_ij(i, j, v::Vector, w::Vector)
    return w[i] * κ_i(i, v, w) * h_i(j, v, w) / (4 * (v[j]^2 - w[i]^2))
end