# Peierls.jl

mutable struct Peierls

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

    function Peierls(x...)
        new(x...)
    end
end

peierls(; kwargs...) = peierls(1; kwargs...)

peierls(α; kwargs...) = peierls(α, 1; kwargs...)

peierls(α, ω; kwargs...) = peierls(α, ω, 1; kwargs...)

peierls(α, ω, J; kwargs...) = peierls(α, ω, J, Inf; kwargs...)

peierls(α, ω, J, β, Ω; kwargs...) = peierls(peierls(α, ω, J, β; kwargs...), Ω; kwargs...)

peierls(h::Peierls; kwargs...) = Peierls((getfield(h, x) for x in fieldnames(Peierls))...)

peierls(h::Array{Peierls}) = Peierls((getfield.(h, x) for x in fieldnames(Peierls))...)

function peierls(α::Number, ω::Number, J::Number, β::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    g = pustrip(peierls_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = peierls_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3) + dims / 2 * log(2π * β) / β, v_guess, w_guess; upper = upper)
    return Peierls(ω * ω0_pu, J * E0_pu, dims, g * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu)
end

function multipeierls(α::AbstractArray, ω::AbstractArray, J::Number, β::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    g = pustrip.(peierls_coupling.(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = sum(peierls_S(v, w, g[j], τ -> β == Inf ? phonon_propagator(τ, ω[j]) : phonon_propagator(τ, ω[j], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims) for j in eachindex(ω))
    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? 1 + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3) + dims / 2 * log(2π * β) / β, v_guess, w_guess; upper = upper)
    return Peierls(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, zero(Float64) .* ω0_pu, zero(Complex) .* ω0_pu)
end

function peierls(α, ω, J, β; verbose = false, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, kwargs...)
    num_α, num_ω, num_J, num_β = length(α), length(ω), length(J), length(β)
    if num_α == num_ω && num_ω > 1 return multipeierls(α, ω, J, β; verbose = verbose, v_guesses = v_guesses, w_guesses = w_guesses, dims = dims, upper = upper) end
    if verbose N, n = num_α * num_ω * num_J * num_β, Threads.Atomic{Int}(1) end
    hstart = peierls(; dims = dims, v_guesses = v_guesses, w_guesses = w_guesses, kwargs...)
    v_guess, w_guess = fill(pustrip.(hstart.v), Threads.nthreads()), fill(pustrip.(hstart.w), Threads.nthreads())
    peierlss = fill(hstart, num_α, num_ω, num_J, num_β)
    Threads.@threads :static for x in CartesianIndices((num_α, num_ω, num_J, num_β))
        id = Threads.threadid()
        if verbose println("\e[KStatics | Threadid: $id | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | α = $(α[x[1]]) [$(x[1])/$num_α] | ω = $(ω[x[2]]) [$(x[2])/$num_ω] | J = $(J[x[3]]) [$(x[3])/$num_J] | β = $(β[x[4]]) [$(x[4])/$num_β]\e[1F"); Threads.atomic_add!(n, 1) end
        @views peierlss[x] = peierls(α[x[1]], ω[x[2]], J[x[3]], β[x[4]]; v_guesses = v_guess[id], w_guesses = w_guess[id], dims = dims, kwargs...)
        v_guess[id], w_guess[id] = pustrip.(peierlss[x].v), pustrip.(peierlss[x].w)
    end
    polaron = peierls(peierlss)
    polaron.α, polaron.g, polaron.d, polaron.ω, polaron.J, polaron.β, polaron.Ω, polaron.Σ = α, reduce_array(polaron.g[:,:,:,1]), dims, pustrip.(ω) * ω0_pu, pustrip.(J) .* E0_pu, pustrip.(β) / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu
    return peierls(polaron)
end

function peierls(h::Peierls, Ω; dims = 3, verbose = false, kwargs...)
    Ω = pustrip.(Ω)
    ω, J, d, g, α, E, v, w, β = [map(y -> pustrip.(y), getfield(h, x)) for x in fieldnames(Peierls)]
    if length(α) == length(ω) return multipeierls(h, Ω; dims = dims, verbose = verbose, kwargs...) end
    num_α, num_ω, num_J, num_β, num_Ω = length(α), length(ω), length(J), length(β), length(Ω)
    if verbose N, n = num_α * num_ω * num_J * num_β * num_Ω, Threads.Atomic{Int}(1) end
    Σ = Array{ComplexF64}(undef, num_α, num_ω, num_J, num_β, num_Ω)
    Threads.@threads for x in CartesianIndices((num_α, num_ω, num_J, num_β, num_Ω))
        if verbose println("\e[KDynamics | Threadid: $(Threads.threadid()) | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | α = $(α[x[1]]) [$(x[1])/$num_α] | ω = $(ω[x[2]]) [$(x[2])/$num_ω] | J = $(J[x[3]]) [$(x[3])/$num_J] | β = $(β[x[4]]) [$(x[4])/$num_β] | Ω = $(Ω[x[5]]) [$(x[5])/$num_Ω]\e[1F"); Threads.atomic_add!(n, 1) end
        @views Σ[x] = peierls_memory(Ω[x[5]], g[x[1],x[2],x[3]], t -> β[x[4]] == Inf ? phonon_propagator(t, ω[x[2]]) : phonon_propagator(t, ω[x[2]], β[x[4]]), t -> β[x[4]] == Inf ? polaron_propagator(t, v[x[1],x[2],x[3],x[4]], w[x[1],x[2],x[3],x[4]]) * J[x[3]] : polaron_propagator(t, v[x[1],x[2],x[3],x[4]], w[x[1],x[2],x[3],x[4]], β[x[4]]) * J[x[3]]; dims = dims) * J[x[3]]
        h.Ω, h.Σ = Ω .* ω0_pu, reduce_array(Σ) .* ω0_pu
    end
    return peierls(h)
end

function multipeierls(α, ω, J, β; verbose = false, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, kwargs...)
    num_J, num_β = length(J), length(β)
    if verbose N, n = num_J * num_β, Threads.Atomic{Int}(1) end
    hstart = peierls(; dims = dims, v_guesses = v_guesses, w_guesses = w_guesses, kwargs...)
    v_guess, w_guess = fill(pustrip.(hstart.v), Threads.nthreads()), fill(pustrip.(hstart.w), Threads.nthreads())
    peierlss = fill(hstart, num_J, num_β)
    Threads.@threads :static for x in CartesianIndices((num_J, num_β))
        id = Threads.threadid()
        if verbose println("\e[KStatics | Threadid: $id | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | J = $(J[x[1]]) [$(x[1])/$num_J] | β = $(β[x[2]]) [$(x[2])/$num_β]\e[1F"); Threads.atomic_add!(n, 1) end
        @views peierlss[x] = peierls(α, ω, J[x[1]], β[x[2]]; v_guesses = v_guess[id], w_guesses = w_guess[id], dims = dims, upper = upper, kwargs...)
        v_guess[id], w_guess[id] = pustrip.(peierlss[x].v), pustrip.(peierlss[x].w)
    end
    polaron = peierls(peierlss)
    polaron.α, polaron.g, polaron.d, polaron.ω, polaron.J, polaron.β, polaron.Ω, polaron.Σ = α, reduce_array(polaron.g[:,1]), dims, pustrip.(ω) * ω0_pu, pustrip.(J) .* E0_pu, pustrip.(β) / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu
    return peierls(polaron)
end

function multipeierls(h::Peierls, Ω; dims = 3, verbose = false, kwargs...)
    Ω = pustrip.(Ω)
    ω, J, d, g, α, E, v, w, β = [map(y -> pustrip.(y), getfield(h, x)) for x in fieldnames(Peierls)]
    num_J, num_β, num_Ω = length(J), length(β), length(Ω)
    if verbose N, n = num_J * num_β * num_Ω, Threads.Atomic{Int}(1) end
    Σ = Array{ComplexF64}(undef, num_J, num_β, num_Ω)
    Threads.@threads for x in CartesianIndices((num_J, num_β, num_Ω))
        if verbose println("\e[KDynamics | Threadid: $(Threads.threadid()) | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | J = $(J[x[1]]) [$(x[1])/$num_J] | β = $(β[x[2]]) [$(x[2])/$num_β] | Ω = $(Ω[x[3]]) [$(x[3])/$num_Ω]\e[1F"); Threads.atomic_add!(n, 1) end
        @views Σ[x] = sum(peierls_memory(Ω[x[3]], g[j], t -> β[x[2]] == Inf ? phonon_propagator(t, ω[j]) : phonon_propagator(t, ω[j], β[x[2]]), t -> β[x[2]] == Inf ? polaron_propagator(t, v[x[1],x[2]], w[x[1],x[2]]) * J[x[1]] : polaron_propagator(t, v[x[1],x[2]], w[x[1],x[2]], β[x[2]]) * J[x[1]]; dims = dims) * J[x[1]] for j in eachindex(ω))
        h.Ω, h.Σ = Ω .* ω0_pu, reduce_array(Σ) .* ω0_pu
    end
    return peierls(h)
end

function peierls(material::Material; kwargs...)
    α = peierls_alpha.(material.g, material.ω_LO, material.J, material.d)
    return peierls(α,puconvert.(material.ω_LO), puconvert.(material.J); dims = material.d, kwargs...)
end

function peierls(material::Material, T; kwargs...)
    α = peierls_alpha.(material.g, material.ω_LO, material.J, material.d)
    return peierls(α, puconvert.(material.ω_LO), puconvert.(material.J), pustrip.(1 ./ (kB_pu .* T)) ./ E0_pu; dims = material.d, kwargs...)
end

function peierls(material::Material, T, Ω; kwargs...)
    α = peierls_alpha.(material.g, material.ω_LO, material.J, material.d)
    return peierls(α, puconvert.(material.ω_LO), puconvert.(material.J), pustrip.(1 ./ (kB_pu .* T)) ./ E0_pu, puconvert.(Ω); dims = material.d, kwargs...)
end

function peierls_alpha(coupling, frequency, transfer_integral, dims)
    
    # Add units
    ω = puconvert(frequency)
    J = puconvert(transfer_integral)
    phonon_energy = ħ_pu * ω

    α = norm(coupling)^2 / (2 * dims * phonon_energy * J)

    return α |> Unitful.NoUnits
end

function peierls_coupling(α, frequency, transfer_integral, dims)

    # Add units
    ω = puconvert(frequency)
    J = puconvert(transfer_integral)
    phonon_energy = ħ_pu * ω

    g = sqrt(2 * dims * α * phonon_energy * J)

    return g |> E0_pu
end

function peierls_S(v, w, coupling, phonon_propagator, polaron_propagator; limits = [0, Inf], dims = 3)
    integral, _ = quadgk(τ -> phonon_propagator(τ) * 1/2 * sqrt(π / polaron_propagator(τ, v, w)) * erf(π * sqrt(polaron_propagator(τ, v, w))) * (1 - exp(-1.0 / polaron_propagator(τ, v, w))), limits...)
    return norm(coupling)^2 * integral / (2π) 
end

function peierls_memory(Ω, coupling, phonon_propagator, polaron_propagator; dims = 3)
    G(t) = polaron_propagator(im * t)
    integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(phonon_propagator(im * t) * (-π/2/G(t) * (1 - exp(-1.0/G(t))) * exp(-π^2 * G(t)) + √π/4/G(t)^(3/2) * (1 - exp(-1.0/G(t))) * erf(π * sqrt(G(t))) + √π / 2 * erf(π * sqrt(G(t))) * exp(-1.0/G(t)) / G(t)^(5/2))), eps(), Inf; rtol=1e-4)
    return norm(coupling)^2 * integral * dims / (2π)^dims
end

function save_peierls(data::Peierls, prefix)

    println("Saving peierls data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "ω", pustrip.(data.ω),
        "J", pustrip.(data.J),
        "d", pustrip.(data.d),
        "g", pustrip.(data.g),
        "α", pustrip.(data.α),
        "E", pustrip.(data.E),
        "v", map(x -> pustrip.(x), data.v),
        "w", map(x -> pustrip.(x), data.w),
        "β", pustrip.(data.β),
        "Ω", pustrip.(data.Ω),
        "Σ", pustrip.(data.Σ)
        )

    println("... peierls data saved.")
end

function load_peierls(peierls_file_path)

    println("Loading peierls data from $peierls_file_path ...")

    data = JLD.load("$peierls_file_path")

    peierls = Peierls(
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
    println("... peierls data loaded.")

    return peierls
end
