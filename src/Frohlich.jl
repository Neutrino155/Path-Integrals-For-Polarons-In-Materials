# Frohlich.jl

mutable struct Frohlich

    # Hamiltonian parameters
    ω   # Phonon frequency
    Mₖ  # Coupling energy
    α   # Dimensionless coupling

    # Polaron properties
    E   # Free Energy
    v   # Variational parameter
    w   # Variational parameter
    β   # Thermodynamic temperature
    Ω   # Electric field frequency
    Σ   # Memory function 

    function Frohlich(x...)
        new(x...)
    end
end

frohlich(; kwargs...) = frohlich(1; kwargs...)

frohlich(α; kwargs...) = frohlich(α, 1; kwargs...)

frohlich(α, ω; kwargs...) = frohlich(α, ω, Inf; kwargs...)

frohlich(α, ω, β, Ω; kwargs...) = frohlich(frohlich(α, ω, β; reduce = false, kwargs...), Ω; kwargs...)

frohlich(f::Frohlich; reduce = true, kwargs...) = reduce ? Frohlich((reduce_array(getfield(f, x)) for x in fieldnames(Frohlich))...) : Frohlich((getfield(f, x) for x in fieldnames(Frohlich))...)

frohlich(f::Array{Frohlich}) = Frohlich((getfield.(f, x) for x in fieldnames(Frohlich))...)

function frohlich(α::Number, ω::Number, β::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), reduce = true, verbose = false)
    ω = pustrip(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = frohlich_S(v, w, Mₖ, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, sqrt(β / 2)])
    v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
    w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
    return Frohlich(ω * ω0_pu, Mₖ * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), reduce = true, verbose = false)
    if length(α) != length(ω) return frohlich(α, ω, β; verbose = verbose, reduce = reduce, v_guesses = v_guesses, w_guesses = w_guesses, mb = mb, dims = dims, upper = upper, lower = lower) end
    ω = pustrip.(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = pustrip.(frohlich_coupling.(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = sum(frohlich_S(v, w, Mₖ[j], τ -> β == Inf ? phonon_propagator(τ, ω[j]) : phonon_propagator(τ, ω[j], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω[j] : polaron_propagator(τ, v, w, β) * ω[j]; dims = dims, limits = [0, sqrt(β / 2)]) for j in eachindex(ω))
    v_guess = v_guesses == false ? sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
    w_guess = w_guesses == false ? 2 + tanh((6 - sum(α)) / 3) + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu)
end

function frohlich(α, ω, β; verbose = false, reduce = true, v_guesses = false, w_guesses = false, kwargs...)
    num_α, num_ω, num_β = length(α), length(ω), length(β)
    if verbose N, n = num_α * num_ω * num_β, Threads.Atomic{Int}(1) end
    fstart = frohlich(; v_guesses = v_guesses, w_guesses = w_guesses, kwargs...)
    v_guess, w_guess = fill(pustrip.(fstart.v), Threads.nthreads()), fill(pustrip.(fstart.w), Threads.nthreads())
    frohlichs = fill(fstart, num_α, num_ω, num_β)
    Threads.@threads :static for x in CartesianIndices((num_α, num_ω, num_β))
        id = Threads.threadid()
        if verbose println("\e[KStatics | Threadid: $id | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | α = $(α[x[1]]) [$(x[1])/$num_α] | ω = $(ω[x[2]]) [$(x[2])/$num_ω] | β = $(β[x[3]]) [$(x[3])/$num_β]\e[1F"); Threads.atomic_add!(n, 1) end
        @views frohlichs[x] = frohlich(α[x[1]], ω[x[2]], β[x[3]]; v_guesses = v_guess[id], w_guesses = w_guess[id], kwargs...)
        v_guess[id], w_guess[id] = pustrip.(frohlichs[x].v), pustrip.(frohlichs[x].w)
    end
    polaron = frohlich(frohlichs)
    polaron.α, polaron.Mₖ, polaron.ω, polaron.β, polaron.Ω, polaron.Σ = α, reduce_array(polaron.Mₖ), pustrip.(ω) * ω0_pu, pustrip.(β) / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu
    return frohlich(polaron; reduce = reduce)
end

function frohlich(f::Frohlich, Ω; dims = 3, verbose = false, kwargs...)
    Ω = pustrip.(Ω)
    ω, Mₖ, α, E, v, w, β = [map(y -> pustrip.(y), getfield(f, x)) for x in fieldnames(Frohlich)]
    if length(α) == length(ω) return multifrohlich(f, Ω; dims = dims, verbose = verbose, kwargs...) end
    num_α, num_ω, num_β, num_Ω = length(α), length(ω), length(β), length(Ω)
    if verbose N, n = num_α * num_ω * num_β * num_Ω, Threads.Atomic{Int}(1) end
    Σ = Array{ComplexF64}(undef, num_α, num_ω, num_β, num_Ω)
    Threads.@threads :static for x in CartesianIndices((num_α, num_ω, num_β, num_Ω))
        if verbose println("\e[KDynamics | Threadid: $(Threads.threadid()) | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | α = $(α[x[1]]) [$(x[1])/$num_α] | ω = $(ω[x[2]]) [$(x[2])/$num_ω] | β = $(β[x[3]]) [$(x[3])/$num_β] | Ω = $(Ω[x[4]]) [$(x[4])/$num_Ω]\e[1F"); Threads.atomic_add!(n, 1) end
        @views Σ[x] = frohlich_memory(Ω[x[4]], Mₖ[x[1],x[2]], t -> β[x[3]] == Inf ? phonon_propagator(t, ω[x[2]]) : phonon_propagator(t, ω[x[2]], β[x[3]]), t -> β[x[3]] == Inf ? polaron_propagator(t, v[x[1],x[2],x[3]], w[x[1],x[2],x[3]]) * ω[x[2]] : polaron_propagator(t, v[x[1],x[2],x[3]], w[x[1],x[2],x[3]], β[x[3]]) * ω[x[2]]; dims = dims) * ω[x[2]]
        f.Ω, f.Σ = Ω .* ω0_pu, reduce_array(Σ) .* ω0_pu
    end
    return frohlich(f)
end

function multifrohlich(f::Frohlich, Ω; dims = 3, verbose = false, kwargs...)
    Ω = pustrip.(Ω)
    ω, Mₖ, α, E, v, w, β = [map(y -> pustrip.(y), getfield(f, x)) for x in fieldnames(Frohlich)]
    num_β, num_Ω = length(β), length(Ω)
    if verbose N, n = num_β * num_Ω, Threads.Atomic{Int}(1) end
    Σ = Array{ComplexF64}(undef, num_β, num_Ω)
    Threads.@threads :static for x in CartesianIndices((num_β, num_Ω))
        if verbose println("\e[KDynamics | Threadid: $(Threads.threadid()) | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | β = $(β[x[1]]) [$(x[1])/$num_β] | Ω = $(Ω[x[2]]) [$(x[2])/$num_Ω]\e[1F"); Threads.atomic_add!(n, 1) end
        @views Σ[x] = sum(frohlich_memory(Ω[x[2]], Mₖ[j], t -> β[x[1]] == Inf ? phonon_propagator(t, ω[j]) : phonon_propagator(t, ω[j], β[x[1]]), t -> β[x[1]] == Inf ? polaron_propagator(t, v[x[1]], w[x[1]]) * ω[j] : polaron_propagator(t, v[x[1]], w[x[1]], β[x[1]]) * ω[j]; dims = dims) * ω[j] for j in eachindex(ω))
        f.Ω, f.Σ = Ω .* ω0_pu, reduce_array(Σ) .* ω0_pu
    end
    return frohlich(f)
end

function frohlich(material::Material; kwargs...)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO; mb = material.mb, kwargs...)
end

function frohlich(material::Material, T; kwargs...)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO, pustrip.(1 ./ (kB_pu .* T)); mb = material.mb, kwargs...)
end

function frohlich(material::Material, T, Ω; kwargs...)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO, pustrip.(1 ./ (kB_pu .* T)), pustrip.(Ω); mb = material.mb, kwargs...)
end
function frohlich_alpha(optical_dielectric, static_dielectic, frequency, effective_mass)
    
    # Add units
    ω = puconvert(frequency)
    mb = puconvert(effective_mass)
    polaron_radius = sqrt(ħ_pu / (2 * mb * ω))
    phonon_energy = ħ_pu * ω
    pekar_factor = (1/puconvert(optical_dielectric) - 1/puconvert(static_dielectic)) / (4π)

    α = (1/2) * e_pu^2 * pekar_factor / (phonon_energy * polaron_radius)

    return α |> Unitful.NoUnits
end

function frohlich_alpha(optical_dielectric, static_dielectic, ionic_dielectric, frequency, effective_mass)
    
    # Add units
    ω = puconvert(frequency)
    mb = puconvert(effective_mass)
    polaron_radius = sqrt(ħ_pu / (2 * mb * ω))
    phonon_energy = ħ_pu * ω
    pekar_factor = puconvert(ionic_dielectric) / (4π * puconvert(static_dielectic) * puconvert(optical_dielectric))

    α = (1/2) * e_pu^2 * pekar_factor / (phonon_energy * polaron_radius)

    return α |> Unitful.NoUnits
end

function ϵ_ionic_mode(frequency, ir_activity, volume)

    # Add units
    ω_j = puconvert(frequency)
    ir_activity = puconvert(ir_activity)
    volume = puconvert(volume)

    # Dielectric contribution from a single ionic phonon mode
    ϵ_mode = ir_activity / (3 * volume * ω_j^2)

    return ϵ_mode |> punit(Unitful.ϵ0)
end

function frohlich_coupling(q, α, frequency, effective_mass; dims = 3)

    # Add units
    ω = puconvert(frequency)
    mb = puconvert(effective_mass)
    phonon_energy = ħ_pu * ω
    q = puconvert(q)
    polaron_radius = sqrt(ħ_pu / (2 * 1m0_pu * 1ω0_pu)) / r0_pu

    Mₖ = im * phonon_energy / q^((dims - 1) / 2) * sqrt(α * polaron_radius * gamma((dims - 1) / 2) * (2√π)^(dims - 1)) 

    return Mₖ |> E0_pu
end

frohlich_coupling(α, frequency, effective_mass; dims = 3) = frohlich_coupling(1, α, frequency, effective_mass; dims = dims)

function frohlich_S(v, w, coupling, phonon_propagator, polaron_propagator; limits = [0, Inf], dims = 3)
    integral, _ = quadgk(τ -> 2 * τ * phonon_propagator(τ^2) / sqrt(polaron_propagator(τ^2, v, w)), limits...)
    return norm(coupling)^2 * ball_surface(dims) / (2π)^dims * sqrt(π / 2) * integral
end

function frohlich_memory(Ω, coupling, phonon_propagator, polaron_propagator; dims = 3)
    integral, _ = quadgk(t -> 2 * t * (1 - exp(im * Ω * t^2)) / Ω * imag(phonon_propagator(im * t^2) / polaron_propagator(im * t^2)^(3/2)), 0, Inf, rtol=1e-4)
    return 2 / dims * norm(coupling)^2 * ball_surface(dims) / (2π)^dims * sqrt(π / 2) * integral
end

function save_frohlich(data::Frohlich, prefix)

    println("Saving Frohlich data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "ω", pustrip.(data.ω),
        "Mₖ", pustrip.(data.Mₖ),
        "α", pustrip.(data.α),
        "E", pustrip.(data.E),
        "v", pustrip.(data.v),
        "w", pustrip.(data.w),
        "β", pustrip.(data.β),
        "Ω", pustrip.(data.Ω),
        "Σ", pustrip.(data.Σ)
        )

    println("... Frohlich data saved.")
end

function load_frohlich(frohlich_file_path)

    println("Loading Frohlich data from $frohlich_file_path ...")

    data = JLD.load("$frohlich_file_path")

    frohlich = Frohlich(
        data["ω"] .* ω0_pu,
        data["Mₖ"] .* E0_pu,
        data["α"],
        data["E"] .* E0_pu,
        data["v"] .* ω0_pu,
        data["w"] .* ω0_pu,
        data["β"] ./ E0_pu,
        data["Ω"] .* ω0_pu,
        data["Σ"] .* ω0_pu
    )
    println("... Frohlich data loaded.")

    return frohlich
end
