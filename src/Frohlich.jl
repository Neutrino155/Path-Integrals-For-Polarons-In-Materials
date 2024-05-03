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
        new(reduce_array.(x)...)
    end
end

frohlich(; kwargs...) = frohlich(1; kwargs...)

frohlich(α; kwargs...) = frohlich(α, 1; kwargs...)

frohlich(α, ω; kwargs...) = frohlich(α, ω, Inf; kwargs...)

frohlich(α, ω, β; kwargs...) = frohlich(α, ω, β, eps(); kwargs...)

function frohlich(α::Number, ω::Number, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = [Inf, Inf], verbose = false)
    ω = pustrip(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = frohlich_S(v, w, Mₖ, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, β / 2])
    v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
    w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
    Σ = frohlich_memory(Ω, Mₖ, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * ω : polaron_propagator(t, v, w, β) * ω; dims = dims) * ω
    return Frohlich(ω * ω0_pu, Mₖ * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, Ω * ω0_pu, Σ * ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = [Inf, Inf], verbose = false)
    ω = pustrip.(ω)
    β = pustrip.(β * ħ_pu / E0_pu) 
    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))
    S(v, w) = sum(frohlich_S(v, w, Mₖ[i], τ -> β == Inf ? phonon_propagator(τ, ω[i]) : phonon_propagator(τ, ω[i], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω[i] : polaron_propagator(τ, v, w, β) * ω[i]; dims = dims, limits = [0, β / 2]) for i in eachindex(ω))
    v_guess = v_guesses == false ? sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
    w_guess = w_guesses == false ? 2 + tanh((6 - sum(α)) / 3) + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
    Σ = sum(frohlich_memory(Ω, Mₖ[i], t -> β == Inf ? phonon_propagator(t, ω[i]) : phonon_propagator(t, ω[i], β), t -> β == Inf ? polaron_propagator(t, v, w) * ω[i] : polaron_propagator(t, v, w, β) * ω[i]; dims = dims) * ω[i] for i in eachindex(ω))
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α, ω, β, Ω; v_guesses = false, w_guesses = false, verbose = false, upper = [Inf, Inf], kwargs...)
    num_α = length(α)
    num_β = length(β)
    num_Ω = length(Ω)
    if verbose N, n = num_α * num_β * num_Ω, 1 end
    Mₖ = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Array{Any}(undef, num_α, num_β, num_Ω)
    @inbounds for i in 1:num_α, j in 1:num_β, k in 1:num_Ω
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
        f = frohlich(α[i], ω, β[j], Ω[k]; v_guesses = v_guesses, w_guesses = w_guesses, upper = upper, kwargs...)
        Mₖ[i], E[i,j], v[i,j], w[i,j], Σ[i,j,k] = f.Mₖ, f.E, f.v, f.w, f.Σ
        v_guesses, w_guesses = pustrip.(v[i,j]), pustrip.(w[i,j])
        upper = pustrip.([maximum(v[i,j]) * 4, maximum(v[i,j]) * 4])
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω, Mₖ, α, E, v, w, pustrip.(β), pustrip.(Ω), Σ)
end

function frohlich(material::Material; kwargs...)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO; kwargs...)
end

function frohlich(material::Material, T; kwargs...)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO, pustrip.(1 ./ (kB_pu .* T)); kwargs...)
end

function frohlich(material::Material, T, Ω; kwargs...)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO, pustrip.(1 ./ (kB_pu .* T)), pustrip.(Ω); kwargs...)
end

function optimal_frohlich(αωβΩ...; rtol = 1e-5, verbose = false, kwargs...) 
    if verbose println("# Fictitious Particles = 1") end
    fs = []
    f = frohlich(αωβΩ...; verbose = true, kwargs...)
    append!(fs, [f])
    energy = f.E
    err = 1
    if verbose n = 2 end
    while all(err .> rtol)
        if verbose println("# Fictitious Particles = $n") end
        v_guesses = length(f.v) > 1 ? pustrip.(vcat(f.v[1], maximum(f.v[1]) * 2)) : pustrip.(vcat(f.v, maximum(f.v) * 2))
        w_guesses = length(f.w) > 1 ? pustrip.(vcat(f.w[1], maximum(f.w[1]) * 2)) : pustrip.(vcat(f.w, maximum(f.w) * 2))
        upper = pustrip.([maximum(v_guesses) * 4, maximum(v_guesses) * 4])
        f = frohlich(αωβΩ...; v_guesses = v_guesses, w_guesses = w_guesses, upper = upper, verbose = verbose, kwargs...)
        append!(fs, f)
        new_energy = f.E
        err = (new_energy .- energy) ./ energy
        energy = new_energy
        if verbose n += 1 end
    end
    return fs
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
    polaron_radius = sqrt(ħ_pu / (2 * mb * ω)) / r0_pu

    Mₖ = im * phonon_energy / q^((dims - 1) / 2) * sqrt(α * polaron_radius * gamma((dims - 1) / 2) * (2√π)^(dims - 1))

    return Mₖ |> E0_pu
end

frohlich_coupling(α, frequency, effective_mass; dims = 3) = frohlich_coupling(1, α, frequency, effective_mass; dims = dims)

function frohlich_S(v, w, coupling, phonon_propagator, polaron_propagator; limits = [0, Inf], dims = 3)
    integral, _ = quadgk(τ -> phonon_propagator(τ) / sqrt(polaron_propagator(τ, v, w)), limits...)
    return norm(coupling)^2 * ball_surface(dims) / (2π)^dims * sqrt(π / 2) * integral
end

function frohlich_memory(Ω, coupling, phonon_propagator, polaron_propagator; dims = 3)
    integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(phonon_propagator(im * t) / polaron_propagator(im * t)^(3/2)), 0, Inf)
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
