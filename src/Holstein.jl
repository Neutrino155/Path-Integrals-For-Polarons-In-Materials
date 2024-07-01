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
        new(x...)
    end
end

holstein(; kwargs...) = holstein(1; kwargs...)

holstein(α; kwargs...) = holstein(α, 1; kwargs...)

holstein(α, ω; kwargs...) = holstein(α, ω, 1; kwargs...)

holstein(α, ω, J; kwargs...) = holstein(α, ω, J, Inf; kwargs...)

holstein(α, ω, J, β, Ω; kwargs...) = holstein(holstein(α, ω, J, β; kwargs...), Ω; kwargs...)

holstein(h::Holstein; kwargs...) = Holstein((getfield(h, x) for x in fieldnames(Holstein))...)

holstein(h::Array{Holstein}) = Holstein((getfield.(h, x) for x in fieldnames(Holstein))...)

function holstein(α::Number, ω::Number, J::Number, β::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = holstein_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::Number, β::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    if length(α) != length(ω) return holstein(α, ω, J, β; verbose = verbose, reduce = reduce, v_guesses = v_guesses, w_guesses = w_guesses, dims = dims, upper = upper) end
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    g = pustrip.(holstein_coupling.(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = sum(holstein_S(v, w, g[j], τ -> β == Inf ? phonon_propagator(τ, ω[j]) : phonon_propagator(τ, ω[j], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims) for j in eachindex(ω))
    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? 1 + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu)
end

function holstein(α, ω, J, β; verbose = false, dims = 3, v_guesses = false, w_guesses = false, kwargs...)
    num_α, num_ω, num_J, num_β = length(α), length(ω), length(J), length(β)
    if verbose N, n = num_α * num_ω * num_J * num_β, Threads.Atomic{Int}(1) end
    hstart = holstein(; dims = dims, v_guesses = v_guesses, w_guesses = w_guesses, kwargs...)
    v_guess, w_guess = fill(pustrip.(hstart.v), Threads.nthreads()), fill(pustrip.(hstart.w), Threads.nthreads())
    holsteins = fill(hstart, num_α, num_ω, num_J, num_β)
    Threads.@threads :static for ijkl in CartesianIndices((num_α, num_ω, num_J, num_β))
        id = Threads.threadid()
        if verbose println("\e[KStatics | Threadid: $id | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | α = $(α[ijkl[1]]) [$(ijkl[1])/$num_α] | ω = $(ω[ijkl[2]]) [$(ijkl[2])/$num_ω] | J = $(J[ijkl[3]]) [$(ijkl[3])/$num_J] | β = $(β[ijkl[4]])\e[1F"); Threads.atomic_add!(n, 1) end
        @views holsteins[ijkl] = holstein(α[ijkl[1]], ω[ijkl[2]], J[ijkl[3]], β[ijkl[4]]; v_guesses = v_guess[id], w_guesses = w_guess[id], dims = dims, kwargs...)
        v_guess[id], w_guess[id] = pustrip.(holsteins[ijkl].v), pustrip.(holsteins[ijkl].w)
    end
    polaron = holstein(holsteins)
    polaron.α, polaron.g, polaron.d, polaron.ω, polaron.J, polaron.β, polaron.Ω, polaron.Σ = α, reduce_array(polaron.g), dims, pustrip.(ω) * ω0_pu, pustrip.(J) .* E0_pu, pustrip.(β) / E0_pu, zero(Float64) * ω0_pu, zero(Complex) * ω0_pu
    return holstein(polaron)
end

function holstein(h::Holstein, Ω; dims = 3, verbose = false, kwargs...)
    Ω = pustrip.(Ω)
    ω, J, d, g, α, E, v, w, β = [map(y -> pustrip.(y), getfield(h, x)) for x in fieldnames(Holstein)]
    if length(α) == length(ω) return multiholstein(h, Ω; dims = dims, verbose = verbose, kwargs...) end
    num_α, num_ω, num_J, num_β, num_Ω = length(α), length(ω), length(J), length(β), length(Ω)
    if verbose N, n = num_α * num_ω * num_J * num_β * num_Ω, Threads.Atomic{Int}(1) end
    Σ = Array{ComplexF64}(undef, num_α, num_ω, num_J, num_β, num_Ω)
    Threads.@threads for ijklm in CartesianIndices((num_α, num_ω, num_J, num_β, num_Ω))
        if verbose println("\e[KDynamics | Threadid: $(Threads.threadid()) | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | α = $(α[ijklm[1]]) [$(ijklm[1])/$num_α] | ω = $(ω[ijklm[2]]) [$(ijklm[2])/$num_ω] | J = $(J[ijklm[3]]) [$(ijklm[3])/$num_J] | β = $(β[ijklm[4]]) [$(ijklm[4])/$num_β] | Ω = $(Ω[ijklm[5]]) [$(ijklm[5])/$num_Ω]\e[1F"); Threads.atomic_add!(n, 1) end
        @views Σ[ijklm] = holstein_memory(Ω[ijklm[5]], g[ijklm[1],ijklm[2],ijklm[3]], t -> β[ijklm[4]] == Inf ? phonon_propagator(t, ω[ijklm[2]]) : phonon_propagator(t, ω[ijklm[2]], β[ijklm[4]]), t -> β[ijklm[4]] == Inf ? polaron_propagator(t, v[ijklm[1],ijklm[2],ijklm[3],ijklm[4]], w[ijklm[1],ijklm[2],ijklm[3],ijklm[4]]) * J[ijklm[3]] : polaron_propagator(t, v[ijklm[1],ijklm[2],ijklm[3],ijklm[4]], w[ijklm[1],ijklm[2],ijklm[3],ijklm[4]], β[ijklm[4]]) * J[ijklm[3]]; dims = dims) * J[ijklm[3]]
        h.Ω, h.Σ = Ω .* ω0_pu, reduce_array(Σ) .* ω0_pu
    end
    return holstein(h)
end

function multiholstein(h::Holstein, Ω; dims = 3, verbose = false, kwargs...)
    Ω = pustrip.(Ω)
    ω, J, d, g, α, E, v, w, β = [map(y -> pustrip.(y), getfield(h, x)) for x in fieldnames(Holstein)]
    num_J, num_β, num_Ω = length(J), length(β), length(Ω)
    if verbose N, n = num_J * num_β * num_Ω, Threads.Atomic{Int}(1) end
    Σ = Array{ComplexF64}(undef, num_J, num_β, num_Ω)
    Threads.@threads for x in CartesianIndices((num_J, num_β, num_Ω))
        if verbose println("\e[KDynamics | Threadid: $(Threads.threadid()) | $(n[])/$N ($(round(n[]/N*100, digits=1)) %)] | J = $(J[x[1]]) [$(x[1])/$num_J] | β = $(β[x[2]]) [$(x[2])/$num_β] | Ω = $(Ω[x[3]]) [$(x[3])/$num_Ω]\e[1F"); Threads.atomic_add!(n, 1) end
        @views Σ[x] = sum(holstein_memory(Ω[x[3]], g[j][x[1]], t -> β[x[2]] == Inf ? phonon_propagator(t, ω[j]) : phonon_propagator(t, ω[j], β[x[2]]), t -> β[x[2]] == Inf ? polaron_propagator(t, v[x[1],x[2]], w[x[1],x[2]]) * J[x[1]] : polaron_propagator(t, v[x[1],x[2]], w[x[1],x[2]], β[x[2]]) * J[x[1]]; dims = dims) * J[x[1]] for j in eachindex(ω))
        h.Ω, h.Σ = Ω .* ω0_pu, reduce_array(Σ) .* ω0_pu
    end
    return holstein(h)
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
    integral, _ = quadgk(t -> (1 - exp(im * Ω * t)) / Ω * imag(phonon_propagator(im * t) * (√π / 2 * erf(π * sqrt(polaron_propagator(im * t))) / polaron_propagator(im * t)^(3/2) - π * exp(-π^2 * polaron_propagator(im * t)) / polaron_propagator(im * t)) * (√π * erf(π * sqrt(polaron_propagator(im * t))) / sqrt(polaron_propagator(im * t)))^(dims - 1)), 0, Inf, rtol=1e-4)
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
        "v", map(x -> pustrip.(x), data.v),
        "w", map(x -> pustrip.(x), data.w),
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
