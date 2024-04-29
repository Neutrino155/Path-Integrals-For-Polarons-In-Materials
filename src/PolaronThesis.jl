module PolaronThesis

# Import packages
using LinearAlgebra, SpecialFunctions, QuadGK, Optim, Unitful, JLD
using Unitful: @unit, Dimension, Dimensions, NoDims, NoUnits, Units, dimension, uconvert, ustrip

# Include files
include("Units.jl")
include("Material.jl")
include("Feynman.jl")
include("Frohlich.jl")
include("Holstein.jl")

# Export functions Å
export puconvert, punit, pustrip, m0_pu, e_pu, ħ_pu, kB_pu, ω0_pu, r0_pu, E0_pu
export Frohlich, frohlich, frohlich_alpha, frohlich_coupling, frohlich_S, frohlich_memory, save_frohlich, load_frohlich
export Holstein, holstein, holstein_alpha, holstein_coupling, holstein_S, holstein_memory, save_holstein, load_holstein
export polaron_propagator, S₀, E₀
export ball_surface, phonon_propagator, variation, reduce_array, ϵ_ionic_mode
export Material, material, FrohlichMaterial, frohlichmaterial, save_frohlich_material, load_frohlich_material, HolsteinMaterial, holsteinmaterial, save_holstein_material, load_holstein_material

# Universal functions

ball_surface(n) = 2 * π^(n/2) / gamma(n/2)

phonon_propagator(τ, ω) = exp(-τ * ω)
phonon_propagator(τ, ω, β) = 1 / (exp(β * ω) - 1) > eps() ? 1 / (exp(β * ω) - 1) * exp(τ * ω) + (1 + 1 / (exp(β * ω) - 1)) * exp(-τ * ω) : exp(-τ * ω)

function variation(energy, initial_v, initial_w; lower = [0, 0], upper = [Inf, Inf], warn = true)

    Δv = initial_v - initial_w # defines a constraint, so that v>w
    initial = [Δv + eps(), initial_w]

    f(x) = energy(x[1] + x[2], x[2])[1]

    solution = Optim.optimize(
        Optim.OnceDifferentiable(f, initial; autodiff=:forward),
        lower,
        upper,
        initial,
        Fminbox(BFGS()),
        Optim.Options(g_tol = eps())
    )

    Δv, w = Optim.minimizer(solution) # v=Δv+w
    energy_minimized = energy(Δv + w, w)

    if !Optim.converged(solution) && warn
        @warn "Failed to converge. v = $(Δv .+ w), w = $w, E = $energy_minimized"
    end

    return Δv + w, w, energy_minimized...
end

variation(energy; lower = [0, 0], upper = [Inf, Inf], warn = true) = variation(energy, 4, 2; lower = lower, upper = upper, warn = warn)

reduce_array(a) = length(a) == 1 ? only(a) : Array(dropdims(a, dims=tuple(findall(size(a) .== 1)...)))

# Register newly defined units with Unitful
Unitful.register(PolaronThesis)
__init__() = Unitful.register(PolaronThesis)

end # module PolaronThesis
