# Holstein.jl

mutable struct Holstein

    # Hamiltonian parameters
    ω
    J
    d
    g
    α

    # Polaron properties
    E
    v
    w
    β
    Ω
    Σ

    function Holstein(x...)
        new(reduce_array.(x)...)
    end
end

holstein(; kwargs...) = hosltein(1; kwargs...)

holstein(α; kwargs...) = holstein(α, 1; kwargs...)

holstein(α, ω; kwargs...) = holstein(α, ω, 1; kwargs...)

holstein(α, ω, J; kwargs...) = holstein(α, ω, J, Inf; kwargs...)

holstein(α, ω, J, β; kwargs...) = holstein(α, ω, J, β, eps(); kwargs...)

function holstein(α::Number, ω::Number, J::Number, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = [Inf, Inf], verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)
    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = holstein_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
    Σ = holstein_memory(Ω, g, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * J : polaron_propagator(t, v, w, β) * J; dims = dims) * J
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α, ω::Number, J::Number, β, Ω; dims = 3, v_guesses = false, w_guesses = false, upper = [Inf, Inf], verbose = false, kwargs...)
    num_α = length(α)
    num_β = length(β)
    num_Ω = length(Ω)
    if verbose N, n = num_α * num_β * num_Ω, 1 end
    g = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Array{Any}(undef, num_α, num_β, num_Ω)
    @inbounds for i in 1:num_α, j in 1:num_β, k in 1:num_Ω
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
        h = holstein(α[i], ω, J, β[j], Ω[k]; v_guesses = v_guesses, w_guesses = w_guesses, upper = upper, dims = dims, kwargs...)
        g[i], v[i,j], w[i,j], E[i,j], Σ[i,j,k] = h.g, h.v, h.w, h.E, h.Σ
        v_guesses, w_guesses = pustrip.(v[i,j]), pustrip.(w[i,j])
        upper = pustrip.([maximum(v[i,j]) * 2 + 2 * ω0_pu, maximum(v[i,j]) * 2 + 2 * ω0_pu])
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω, J, dims, g, α, E, v, w, β, Ω, Σ)
end

function holstein(material::Material; kwargs...)
    α = holstein_alpha.(material.g, material.ω_LO, material.J, material.d)
    return holstein(α,puconvert.(material.ω_LO), puconvert.(material.J); dims = material.d, kwargs...)
end

function holstein(material::Material, T; kwargs...)
    α = holstein_alpha.(material.g, material.ω_LO, material.J, material.d)
    return holstein(α, puconvert.(material.ω_LO), puconvert.(material.J), pustrip.(1 ./ (kB_pu .* T)) ./ E0_pu; dims = material.d, kwargs...)
end

function holstein(material::Material, T, Ω; kwargs...)
    α = holstein_alpha.(material.g, material.ω_LO, material.J, material.d)
    return holstein(α, puconvert.(material.ω_LO), puconvert.(material.J), pustrip.(1 ./ (kB_pu .* T)) ./ E0_pu, puconvert.(Ω); dims = material.d, kwargs...)
end

function optimal_holstein(αωJβΩ...; rtol = 1e-5, verbose = false, kwargs...) 
    if verbose println("# Fictitious Particles = 1") end
    h = holstein(αωJβΩ...; verbose = verbose, kwargs...)
    energy = h.E
    err = 1
    if verbose n = 2 end
    while all(err .> rtol)
        if verbose println("# Fictitious Particles = $n") end
        v_guesses = length(h.v) > 1 ? pustrip.(vcat(h.v[1], maximum(h.v[1]) * 2)) : pustrip.(vcat(h.v, maximum(h.v) * 2))
        w_guesses = length(h.w) > 1 ? pustrip.(vcat(h.w[1], maximum(h.w[1]) * 2)) : pustrip.(vcat(h.w, maximum(h.w) * 2))
        upper = pustrip.([maximum(v_guesses) * 4, maximum(v_guesses) * 4])
        h = holstein(αωJβΩ...; v_guesses = v_guesses, w_guesses = w_guesses, upper = upper, verbose = verbose, kwargs...)
        new_energy = h.E
        err = (new_energy .- energy) ./ energy 
        energy = new_energy
        if verbose n += 1 end
    end
    return h
end
function holstein_alpha(coupling, frequency, transfer_integral, dims)
    
    # Add units
    ω = puconvert(frequency)
    J = puconvert(transfer_integral)
    phonon_energy = ħ_pu * ω

    α = norm(coupling)^2 / (2 * dims * phonon_energy * J)

    return α |> Unitful.NoUnits
end

function holstein_coupling(α, frequency, transfer_integral, dims)

    # Add units
    ω = puconvert(frequency)
    J = puconvert(transfer_integral)
    phonon_energy = ħ_pu * ω

    g = sqrt(2 * dims * α * phonon_energy * J)

    return g |> E0_pu
end

function holstein_S(v, w, coupling, phonon_propagator, polaron_propagator; limits = [0, Inf], dims = 3)
    integral, _ = quadgk(τ -> phonon_propagator(τ) * (sqrt(π / polaron_propagator(τ, v, w)) * erf(π * sqrt(polaron_propagator(τ, v, w))))^dims, limits...)
    return norm(coupling)^2 * integral / (2π)^dims / 2
end

function holstein_memory(Ω, coupling, phonon_propagator, polaron_propagator; dims = 3)
    integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(phonon_propagator(im * t) * (√π / 2 * erf(π * sqrt(polaron_propagator(im * t))) / polaron_propagator(im * t)^(3/2) - π * exp(-π^2 * polaron_propagator(im * t)) / polaron_propagator(im * t)) * (√π * erf(π * sqrt(polaron_propagator(im * t))) / sqrt(polaron_propagator(im * t)))^(dims - 1)), 0, Inf)
    return norm(coupling)^2 * integral * dims / (2π)^dims
end

function save_holstein(data::Holstein, prefix)

    println("Saving holstein data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "ω", pustrip.(data.ω),
        "J", pustrip.(data.J),
        "d", pustrip.(data.d),
        "g", pustrip.(data.g),
        "α", pustrip.(data.α),
        "E", pustrip.(data.E),
        "v", pustrip.(data.v),
        "w", pustrip.(data.w),
        "β", pustrip.(data.β),
        "Ω", pustrip.(data.Ω),
        "Σ", pustrip.(data.Σ)
        )

    println("... holstein data saved.")
end

function load_holstein(holstein_file_path)

    println("Loading holstein data from $holstein_file_path ...")

    data = JLD.load("$holstein_file_path")

    holstein = Holstein(
        data["ω"] .* ω0_pu,
        data["J"] .* E0_pu,
        data["d"],
        data["g"] .* E0_pu,
        data["α"],
        data["E"] .* E0_pu,
        data["v"] .* ω0_pu,
        data["w"] .* ω0_pu,
        data["β"] ./ E0_pu,
        data["Ω"] .* ω0_pu,
        data["Σ"] .* ω0_pu
    )
    println("... holstein data loaded.")

    return holstein
end
