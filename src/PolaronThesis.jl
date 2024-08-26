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
include("Peierls.jl")

# Export functions Å
export puconvert, punit, pustrip, m0_pu, e_pu, ħ_pu, kB_pu, ω0_pu, r0_pu, E0_pu
export Frohlich, frohlich, frohlich_alpha, frohlich_coupling, frohlich_S, frohlich_memory, save_frohlich, load_frohlich
export Holstein, holstein, holstein_alpha, holstein_coupling, holstein_S, holstein_memory, save_holstein, load_holstein, optimal_holstein
export Peierls, peierls, peierls_alpha, peierls_coupling, peierls_S, peierels_memory, save_peierls, load_peierls, optimal_peierls
export polaron_propagator, S₀, E₀
export ball_surface, phonon_propagator, variation, reduce_array, reshape_array, LogRange, ϵ_ionic_mode
export Material, material, FrohlichMaterial, frohlichmaterial, save_frohlich_material, load_frohlich_material, HolsteinMaterial, holsteinmaterial, save_holstein_material, load_holstein_material

# Universal functions

ball_surface(n) = 2 * π^(n/2) / gamma(n/2)

phonon_propagator(τ, ω) = exp(-τ * ω)
phonon_propagator(τ, ω, β) = 1 / (exp(β * ω) - 1) > eps() ? 1 / (exp(β * ω) - 1) * exp(τ * ω) + (1 + 1 / (exp(β * ω) - 1)) * exp(-τ * ω) : exp(-τ * ω)

function variation(energy, initial_v, initial_w; lower = eps(), upper = Inf, warn = false)

    if length(initial_v) != length(initial_w)
        return error("The number of variational parameters v & w must be equal.")
    end
  
    N_params = length(initial_v)

    params = sort(vcat(initial_w, initial_v))

    Δx = params[2:end] .- params[1:end-1]
    initial = vcat(minimum(initial_w), [Δx .+ eps()]...)
    
    f(x) = energy([sum(x[i] for i in 1:(2 * n)) for n in 1:N_params], [sum(x[i] for i in 1:(2 * n - 1)) for n in 1:N_params])[1]
  
    solution = Optim.optimize(
      Optim.OnceDifferentiable(f, initial; autodiff=:forward),
      lower,
      upper,
      initial,
      Fminbox(BFGS()),
      Optim.Options(g_tol = eps())
    )
  
    var_params = cumsum(Optim.minimizer(solution))
  
    v = reduce_array([var_params[2*n] for n in 1:N_params])
    w = reduce_array([var_params[2*n-1] for n in 1:N_params])

    energy_minimized = reduce_array(Optim.minimum(solution))
  
    if !Optim.converged(solution) && warn
          @warn "Failed to converge. v = $v, w = $w, E = $energy_minimized"
    end
  
    return v, w, energy_minimized...
end

variation(energy; lower = [0, 0], upper = [Inf, Inf], warn = true) = variation(energy, 4, 2; lower = lower, upper = upper, warn = warn)

reduce_array(a) = length(a) == 1 ? only(a) : Array(dropdims(a, dims=tuple(findall(size(a) .== 1)...)))
reshape_array(M) = [M[I][k] for I=CartesianIndices(M),k=eachindex(M[1,1])]
LogRange(i, f, n) = exp.(LinRange(log(i), log(f), n))

# Register newly defined units with Unitful
Unitful.register(PolaronThesis)
__init__() = Unitful.register(PolaronThesis)

end # module PolaronThesis
