# Feynman.jl

function polaron_propagator(τ, v, w)
    c = (v^2 - w^2) / v^3
    result = (1 - v * c) * τ + c * (phonon_propagator(0, v) - phonon_propagator(τ, v))
    return result + sqrt(eps())
end

function polaron_propagator(τ, v, w, β)
    c = (v^2 - w^2) / v^3
    result = c * (phonon_propagator(0, v, β) - phonon_propagator(τ, v, β)) + (1 - c * v) * τ * (1 - τ / β)
    return result + sqrt(eps())
end

S₀(v, w) = -(3 / (4 * v)) * (v^2 - w^2)

S₀(v, w, β) = -3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β)) 

E₀(v, w) = 3 * (v - w) / 2

E₀(v, w, β) = 3 * (v - w) / 2 - 3 / β * (log(v / w) - log(1 - exp(-v * β)) + log(1 - exp(-w * β))) 