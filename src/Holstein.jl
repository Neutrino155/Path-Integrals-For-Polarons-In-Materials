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

holstein(; kwargs...) = holstein(1; kwargs...)

holstein(α; kwargs...) = holstein(α, 1; kwargs...)

holstein(α, ω; kwargs...) = holstein(α, ω, 1; kwargs...)

holstein(α, ω, J; kwargs...) = holstein(α, ω, J, Inf; kwargs...)

holstein(α, ω, J, β; kwargs...) = holstein(α, ω, J, β, eps(); kwargs...)

function holstein(α::Number, ω::Number, J::Number, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)
    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = holstein_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
    Σ = holstein_memory(Ω, g, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * J : polaron_propagator(t, v, w, β) * J; dims = dims) * J
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, Ω * ω0_pu, Σ * ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::Number, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α = length(α)

    if verbose N, n = num_α, 1 end

    g = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Vector{Any}(undef, num_α)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses

    for i in eachindex(α)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α]") end
        g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
        v[i], w[i], E[i] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[i] = holstein_memory(Ω, g[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) * J : polaron_propagator(t, v[i], w[i], β) * J; dims = dims) * J
        v_guess = v_guesses == false ? v[i] : v_guesses
        w_guess = w_guesses == false ? w[i] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::Number, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_ω = length(ω)

    if verbose N, n = num_ω, 1 end

    g = Vector{Any}(undef, num_ω)
    v = Vector{Any}(undef, num_ω)
    w = Vector{Any}(undef, num_ω)
    E = Vector{Any}(undef, num_ω)
    Σ = Vector{Any}(undef, num_ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[p] + 1 / β : w_guesses

    for p in eachindex(ω)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω]") end
        g[p] = pustrip(holstein_coupling(α, ω[p] * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
        v[p], w[p], E[p] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[p] = holstein_memory(Ω, g[p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[p], w[p]) * J : polaron_propagator(t, v[p], w[p], β) * J; dims = dims) * J
        v_guess = v_guesses == false ? v[p] : v_guesses
        w_guess = w_guesses == false ? w[p] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::AbstractArray, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_J = length(J)

    if verbose N, n = num_J, 1 end

    g = Vector{Any}(undef, num_J)
    v = Vector{Any}(undef, num_J)
    w = Vector{Any}(undef, num_J)
    E = Vector{Any}(undef, num_J)
    Σ = Vector{Any}(undef, num_J)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses

    for q in eachindex(J)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | J = $(J[q]) [$q/$num_J]") end
        g[q] = pustrip(holstein_coupling(α, ω * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[q], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[q], w[q], E[q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[q] = holstein_memory(Ω, g[q], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[q], w[q]) * J[q] : polaron_propagator(t, v[q], w[q], β) * J[q]; dims = dims) * J[q]
        v_guess = v_guesses == false ? v[q] : v_guesses
        w_guess = w_guesses == false ? w[q] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω * ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_β = length(β)

    if verbose N, n = num_β, 1 end

    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)

    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β[1] : w_guesses

    for j in eachindex(β)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = holstein_S(v, w, g, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
        v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[j] = holstein_memory(Ω, g, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) * J : polaron_propagator(t, v[j], w[j], β[j]) * J; dims = dims) * J
        v_guess = v_guesses == false ? v[j] : v_guesses
        w_guess = w_guesses == false ? w[j] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g * E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_Ω = length(Ω)

    if verbose N, n = num_Ω, 1 end

    Σ = Vector{Any}(undef, num_Ω)

    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims)) 
    S(v, w) = holstein_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)

    for k in eachindex(Ω)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | Ω = $(Ω[k]) [$k/$num_Ω]") end
        Σ[k] = holstein_memory(Ω[k], g, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * J : polaron_propagator(t, v, w, β) * J; dims = dims) * J
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::Number, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_ω = length(α), length(ω)

    if verbose N, n = num_α * num_ω, 1 end

    g = Matrix{Any}(undef, num_α, num_ω)
    v = Matrix{Any}(undef, num_α, num_ω)
    w = Matrix{Any}(undef, num_α, num_ω)
    E = Matrix{Any}(undef, num_α, num_ω)
    Σ = Matrix{Any}(undef, num_α, num_ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β : w_guesses

    for i in eachindex(α), p in eachindex(ω)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω]") end
        g[i,p] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i,p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
        v[i,p], w[i,p], E[i,p] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[i,p] = holstein_memory(Ω, g[i,p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[i,p], w[i,p]) * J : polaron_propagator(t, v[i,p], w[i,p], β) * J; dims = dims) * J
        v_guess = v_guesses == false ? v[i,p] : v_guesses
        w_guess = w_guesses == false ? w[i,p] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end

    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::AbstractArray, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_J = length(α), length(J)

    if verbose N, n = num_α * num_J, 1 end

    g = Matrix{Any}(undef, num_α, num_J)
    v = Matrix{Any}(undef, num_α, num_J)
    w = Matrix{Any}(undef, num_α, num_J)
    E = Matrix{Any}(undef, num_α, num_J)
    Σ = Matrix{Any}(undef, num_α, num_J)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses

    for i in eachindex(α), q in eachindex(J)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | J = $(J[q]) [$q/$num_J]") end
        g[i,q] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i,q], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[i,q], w[i,q], E[i,q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[i,q] = holstein_memory(Ω, g[i,q], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i,q], w[i,q]) * J[q] : polaron_propagator(t, v[i,q], w[i,q], β) * J[q]; dims = dims) * J[q]
        v_guess = v_guesses == false ? v[i,q] : v_guesses
        w_guess = w_guesses == false ? w[i,q] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω * ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::Number, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_β = length(α), length(β)

    if verbose N, n = num_α * num_β, 1 end

    g = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Matrix{Any}(undef, num_α, num_β)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β[1] : w_guesses

    for i in eachindex(α)
        g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = holstein_S(v, w, g[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
            v[i,j], w[i,j], E[i,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            Σ[i,j] = holstein_memory(Ω, g[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,j], w[i,j]) * J : polaron_propagator(t, v[i,j], w[i,j], β[j]) * J; dims = dims) * J
            v_guess = v_guesses == false ? v[i,j] : v_guesses
            w_guess = w_guesses == false ? w[i,j] : w_guesses
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::Number, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_Ω = length(α), length(Ω)

    if verbose N, n = num_α * num_Ω, 1 end

    g = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Matrix{Any}(undef, num_α, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses

    for i in eachindex(α)
        g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
        v[i], w[i], E[i] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[i] : v_guesses
        w_guess = w_guesses == false ? w[i] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[i,k] = holstein_memory(Ω[k], g[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) * J : polaron_propagator(t, v[i], w[i], β) * J; dims = dims) * J
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::AbstractArray, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_ω, num_J = length(ω), length(J)

    if verbose N, n = num_ω * num_J, 1 end

    g = Matrix{Any}(undef, num_ω, num_J)
    v = Matrix{Any}(undef, num_ω, num_J)
    w = Matrix{Any}(undef, num_ω, num_J)
    E = Matrix{Any}(undef, num_ω, num_J)
    Σ = Matrix{Any}(undef, num_ω, num_J)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β : w_guesses

    for p in eachindex(ω), q in eachindex(J)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J]") end
        g[p,q] = pustrip(holstein_coupling(α, ω[p] * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[p,q], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[p,q], w[p,q], E[p,q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[p,q] = holstein_memory(Ω, g[p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[p,q], w[p,q]) * J[q] : polaron_propagator(t, v[p,q], w[p,q], β) * J[q]; dims = dims) * J[q]
        v_guess = v_guesses == false ? v[p,q] : v_guesses
        w_guess = w_guesses == false ? w[p,q] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::Number, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_ω, num_β = length(ω), length(β)

    if verbose N, n = num_ω * num_β, 1 end

    g = Vector{Any}(undef, num_ω)
    v = Matrix{Any}(undef, num_ω, num_β)
    w = Matrix{Any}(undef, num_ω, num_β)
    E = Matrix{Any}(undef, num_ω, num_β)
    Σ = Matrix{Any}(undef, num_ω, num_β)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for p in eachindex(ω)
        g[p] = pustrip(holstein_coupling(α, ω[p] * ω0_pu, J * E0_pu, dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = holstein_S(v, w, g[p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
            v[p,j], w[p,j], E[p,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            Σ[p,j] = holstein_memory(Ω, g[p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[p,j], w[p,j]) * J : polaron_propagator(t, v[p,j], w[p,j], β[j]) * J; dims = dims) * J
            v_guess = v_guesses == false ? v[p,j] : v_guesses
            w_guess = w_guesses == false ? w[p,j] : w_guesses
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::Number, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_ω, num_Ω = length(ω), length(Ω)

    if verbose N, n = num_ω * num_Ω, 1 end

    g = Vector{Any}(undef, num_ω)
    v = Vector{Any}(undef, num_ω)
    w = Vector{Any}(undef, num_ω)
    E = Vector{Any}(undef, num_ω)
    Σ = Matrix{Any}(undef, num_ω, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β : w_guesses

    for p in eachindex(ω)
        g[p] = pustrip(holstein_coupling(α, ω[p] * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
        v[p], w[p], E[p] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[p] : v_guesses
        w_guess = w_guesses == false ? w[p] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[p,k] = holstein_memory(Ω[k], g[p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[p], w[p]) * J : polaron_propagator(t, v[p], w[p], β) * J; dims = dims) * J
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::AbstractArray, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_J, num_β = length(J), length(β)

    if verbose N, n = num_J * num_β, 1 end

    g = Vector{Any}(undef, num_J)
    v = Matrix{Any}(undef, num_J, num_β)
    w = Matrix{Any}(undef, num_J, num_β)
    E = Matrix{Any}(undef, num_J, num_β)
    Σ = Matrix{Any}(undef, num_J, num_β)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β[1] : w_guesses

    for q in eachindex(J)
        g[q] = pustrip(holstein_coupling(α, ω * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | J = $(J[q]) [$q/$num_J] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = holstein_S(v, w, g[q], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[q,j], w[q,j], E[q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            Σ[q,j] = holstein_memory(Ω, g[q], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[q,j], w[q,j]) * J[q] : polaron_propagator(t, v[q,j], w[q,j], β[j]) * J[q]; dims = dims) * J[q]
            v_guess = v_guesses == false ? v[q,j] : v_guesses
            w_guess = w_guesses == false ? w[q,j] : w_guesses
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω * ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::AbstractArray, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_J, num_Ω = length(J), length(Ω)

    if verbose N, n = num_J * num_Ω, 1 end

    g = Vector{Any}(undef, num_J)
    v = Vector{Any}(undef, num_J)
    w = Vector{Any}(undef, num_J)
    E = Vector{Any}(undef, num_J)
    Σ = Matrix{Any}(undef, num_J, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses

    for q in eachindex(J)
        g[q] = pustrip(holstein_coupling(α, ω * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[q], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[q], w[q], E[q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[q] : v_guesses
        w_guess = w_guesses == false ? w[q] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | J = $(J[q]) [$q/$num_J] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[q,k] = holstein_memory(Ω[k], g[q], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[q], w[q]) * J[q] : polaron_propagator(t, v[q], w[q], β) * J[q]; dims = dims) * J[q]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω * ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_β, num_Ω, = length(β), length(Ω)

    if verbose N, n = num_β * num_Ω, 1 end

    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Matrix{Any}(undef, num_β, num_Ω)

    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β[1] : w_guesses

    for j in eachindex(β)
        S(v, w) = holstein_S(v, w, g, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
        v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[j] : v_guesses
        w_guess = w_guesses == false ? w[j] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[j,k] = holstein_memory(Ω[k], g, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) * J : polaron_propagator(t, v[j], w[j], β[j]) * J; dims = dims) * J
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::AbstractArray, β::Number, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_ω, num_J = length(α), length(ω), length(J)

    if verbose N, n = num_α * num_ω * num_J, 1 end

    g = Array{Any}(undef, num_α, num_ω, num_J)
    v = Array{Any}(undef, num_α, num_ω, num_J)
    w = Array{Any}(undef, num_α, num_ω, num_J)
    E = Array{Any}(undef, num_α, num_ω, num_J)
    Σ = Array{Any}(undef, num_α, num_ω, num_J)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β : w_guesses

    for i in eachindex(α), p in eachindex(ω), q in eachindex(J)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J]") end
        g[i,p,q] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i,p,q], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[i,p,q], w[i,p,q], E[i,p,q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        Σ[i,p,q] = holstein_memory(Ω, g[i,p,q], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[i,p,q], w[i,p,q]) * J[q] : polaron_propagator(t, v[i,p,q], w[i,p,q], β) * J[q]; dims = dims) * J[q]
        v_guess = v_guesses == false ? v[i,p,q] : v_guesses
        w_guess = w_guesses == false ? w[i,p,q] : w_guesses
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::Number, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_ω, num_β = length(α), length(ω), length(β)

    if verbose N, n = num_α * num_ω * num_β, 1 end

    g = Matrix{Any}(undef, num_α, num_ω)
    v = Array{Any}(undef, num_α, num_ω, num_β)
    w = Array{Any}(undef, num_α, num_ω, num_β)
    E = Array{Any}(undef, num_α, num_ω, num_β)
    Σ = Array{Any}(undef, num_α, num_ω, num_β)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for i in eachindex(α), p in eachindex(ω)
        g[i,p] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J * E0_pu, dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = holstein_S(v, w, g[i,p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
            v[i,p,j], w[i,p,j], E[i,p,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            Σ[i,p,j] = holstein_memory(Ω, g[i,p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,p,j], w[i,p,j]) * J : polaron_propagator(t, v[i,p,j], w[i,p,j], β[j]) * J; dims = dims) * J
            v_guess = v_guesses == false ? v[i,p,j] : v_guesses
            w_guess = w_guesses == false ? w[i,p,j] : w_guesses
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::Number, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_ω, num_Ω = length(α), length(ω), length(Ω)

    if verbose N, n = num_α * num_ω * num_Ω, 1 end

    g = Matrix{Any}(undef, num_α, num_ω)
    v = Matrix{Any}(undef, num_α, num_ω)
    w = Matrix{Any}(undef, num_α, num_ω)
    E = Matrix{Any}(undef, num_α, num_ω)
    Σ = Array{Any}(undef, num_α, num_ω, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β : w_guesses

    for i in eachindex(α), p in eachindex(ω)
        g[i,p] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i,p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β) * J; limits = [0, β / 2], dims = dims)
        v[i,p], w[i,p], E[i,p] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[i,p] : v_guesses
        w_guess = w_guesses == false ? w[i,p] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[i,p,k] = holstein_memory(Ω[k], g[i,p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[i,p], w[i,p]) * J : polaron_propagator(t, v[i,p], w[i,p], β) * J; dims = dims) * J
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::AbstractArray, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_J, num_β = length(α), length(J), length(β)

    if verbose N, n = num_α * num_J * num_β, 1 end

    g = Matrix{Any}(undef, num_α, num_J)
    v = Array{Any}(undef, num_α, num_J, num_β)
    w = Array{Any}(undef, num_α, num_J, num_β)
    E = Array{Any}(undef, num_α, num_J, num_β)
    Σ = Array{Any}(undef, num_α, num_J, num_β)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β[1] : w_guesses

    for i in eachindex(α), q in eachindex(J)
        g[i,q] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | J = $(J[q]) [$q/$num_J] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = holstein_S(v, w, g[i,q], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[i,q,j], w[i,q,j], E[i,q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            Σ[i,q,j] = holstein_memory(Ω, g[i,q], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,q,j], w[i,q,j]) * J[q] : polaron_propagator(t, v[i,q,j], w[i,q,j], β[j]) * J[q]; dims = dims) * J[q]
            v_guess = v_guesses == false ? v[i,q,j] : v_guesses
            w_guess = w_guesses == false ? w[i,q,j] : w_guesses
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω * ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::AbstractArray, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_J, num_Ω = length(α), length(J), length(Ω)

    if verbose N, n = num_α * num_J * num_Ω, 1 end

    g = Matrix{Any}(undef, num_α, num_J)
    v = Matrix{Any}(undef, num_α, num_J)
    w = Matrix{Any}(undef, num_α, num_J)
    E = Matrix{Any}(undef, num_α, num_J)
    Σ = Array{Any}(undef, num_α, num_J, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β : w_guesses

    for i in eachindex(α), q in eachindex(J)
        g[i,q] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i,q], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[i,q], w[i,q], E[i,q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[i,q] : v_guesses
        w_guess = w_guesses == false ? w[i,q] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | J = $(J[q]) [$q/$num_J] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[i,q,k] = holstein_memory(Ω[k], g[i,q], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i,q], w[i,q]) * J[q] : polaron_propagator(t, v[i,q], w[i,q], β) * J[q]; dims = dims) * J[q]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω * ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::AbstractArray, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_ω, num_J, num_β = length(ω), length(J), length(β)

    if verbose N, n = num_ω * num_J * num_β, 1 end

    g = Matrix{Any}(undef, num_ω, num_J)
    v = Array{Any}(undef, num_ω, num_J, num_β)
    w = Array{Any}(undef, num_ω, num_J, num_β)
    E = Array{Any}(undef, num_ω, num_J, num_β)
    Σ = Array{Any}(undef, num_ω, num_J, num_β)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for p in eachindex(ω), q in eachindex(J)
        g[p,q] = pustrip(holstein_coupling(α, ω[p] * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = holstein_S(v, w, g[p,q], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[p,q,j], w[p,q,j], E[p,q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            Σ[p,q,j] = holstein_memory(Ω, g[p,q], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[p,q,j], w[p,q,j]) * J[q] : polaron_propagator(t, v[p,q,j], w[p,q,j], β[j]) * J[q]; dims = dims) * J[q]
            v_guess = v_guesses == false ? v[p,q,j] : v_guesses
            w_guess = w_guesses == false ? w[p,q,j] : w_guesses
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::AbstractArray, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_ω, num_J, num_Ω = length(ω), length(J), length(Ω)

    if verbose N, n = num_ω * num_J * num_Ω, 1 end

    g = Matrix{Any}(undef, num_ω, num_J)
    v = Matrix{Any}(undef, num_ω, num_J)
    w = Matrix{Any}(undef, num_ω, num_J)
    E = Matrix{Any}(undef, num_ω, num_J)
    Σ = Array{Any}(undef, num_ω, num_J, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β : w_guesses

    for p in eachindex(ω), q in eachindex(J)
        g[p,q] = pustrip(holstein_coupling(α, ω[p] * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[p,q], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[p,q], w[p,q], E[p,q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[p,q] : v_guesses
        w_guess = w_guesses == false ? w[p,q] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[p,q,k] = holstein_memory(Ω[k], g[p,q], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[p,q], w[p,q]) * J[q] : polaron_propagator(t, v[p,q], w[p,q], β) * J[q]; dims = dims) * J[q]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::Number, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_β, num_Ω, = length(α), length(β), length(Ω)

    if verbose N, n = num_α * num_β * num_Ω, 1 end

    g = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Array{Any}(undef, num_α, num_β, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β[1] : w_guesses

    for i in eachindex(α)
        g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        for j in eachindex(β)
            S(v, w) = holstein_S(v, w, g[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
            v[i,j], w[i,j], E[i,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            v_guess = v_guesses == false ? v[i,j] : v_guesses
            w_guess = w_guesses == false ? w[i,j] : w_guesses
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[i,j,k] = holstein_memory(Ω[k], g[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,j], w[i,j]) * J : polaron_propagator(t, v[i,j], w[i,j], β[j]) * J; dims = dims) * J
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Holstein(ω * ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::Number, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_ω, num_β, num_Ω, = length(ω), length(β), length(Ω)

    if verbose N, n = num_ω * num_β * num_Ω, 1 end

    g = Vector{Any}(undef, num_ω)
    v = Matrix{Any}(undef, num_ω, num_β)
    w = Matrix{Any}(undef, num_ω, num_β)
    E = Matrix{Any}(undef, num_ω, num_β)
    Σ = Array{Any}(undef, num_ω, num_β, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for p in eachindex(ω)
        g[p] = pustrip(holstein_coupling(α, ω[p] * ω0_pu, J * E0_pu, dims))
        for j in eachindex(β)
            S(v, w) = holstein_S(v, w, g[p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
            v_guess = v_guesses == false ? 2 * dims + 1 / β[j] : v_guesses
            w_guess = w_guesses == false ? ω[p] + 1 / β[j] : w_guesses
            v[p,j], w[p,j], E[p,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            v_guess = v_guesses == false ? v[p,j] : v_guesses
            w_guess = w_guesses == false ? w[p,j] : w_guesses
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[p,j,k] = holstein_memory(Ω[k], g[p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[p,j], w[p,j]) * J : polaron_propagator(t, v[p,j], w[p,j], β[j]) * J; dims = dims) * J
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::AbstractArray, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_J, num_β, num_Ω, = length(J), length(β), length(Ω)

    if verbose N, n = num_J * num_β * num_Ω, 1 end

    g = Vector{Any}(undef, num_J)
    v = Matrix{Any}(undef, num_J, num_β)
    w = Matrix{Any}(undef, num_J, num_β)
    E = Matrix{Any}(undef, num_J, num_β)
    Σ = Array{Any}(undef, num_J, num_β, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for q in eachindex(J)
        g[q] = pustrip(holstein_coupling(α, ω * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            S(v, w) = holstein_S(v, w, g[q], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[q,j], w[q,j], E[q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            v_guess = v_guesses == false ? v[q,j] : v_guesses
            w_guess = w_guesses == false ? w[q,j] : w_guesses
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[q,j,k] = holstein_memory(Ω[k], g[q], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[q,j], w[q,j]) * J[q] : polaron_propagator(t, v[q,j], w[q,j], β[j]) * J[q]; dims = dims) * J[q]
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::AbstractArray, β::AbstractArray, Ω::Number; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_ω, num_J, num_β = length(α), length(ω), length(J), length(β)

    if verbose N, n = num_α * num_ω * num_J * num_β, 1 end

    g = Array{Any}(undef, num_α, num_ω, num_J)
    v = Array{Any}(undef, num_α, num_ω, num_J, num_β)
    w = Array{Any}(undef, num_α, num_ω, num_J, num_β)
    E = Array{Any}(undef, num_α, num_ω, num_J, num_β)
    Σ = Array{Any}(undef, num_α, num_ω, num_J, num_β)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for i in eachindex(α), p in eachindex(ω), q in eachindex(J)
        g[i,p,q] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = holstein_S(v, w, g[i,p,q], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[i,p,q,j], w[i,p,q,j], E[i,p,q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            Σ[i,p,q,j] = holstein_memory(Ω, g[i,p,q], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,p,q,j], w[i,p,q,j]) * J[q] : polaron_propagator(t, v[i,p,q,j], w[i,p,q,j], β[j]) * J[q]; dims = dims) * J[q]
            v_guess = v_guesses == false ? v[i,p,q,j] : v_guesses
            w_guess = w_guesses == false ? w[i,p,q,j] : w_guesses
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::AbstractArray, β::Number, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_ω, num_J, num_Ω = length(α), length(ω), length(J), length(Ω)

    if verbose N, n = num_α * num_ω * num_J * num_Ω, 1 end

    g = Array{Any}(undef, num_α, num_ω, num_J)
    v = Array{Any}(undef, num_α, num_ω, num_J)
    w = Array{Any}(undef, num_α, num_ω, num_J)
    E = Array{Any}(undef, num_α, num_ω, num_J)
    Σ = Array{Any}(undef, num_α, num_ω, num_J, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β : w_guesses

    for i in eachindex(α), p in eachindex(ω), q in eachindex(J)
        g[i,p,q] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J[q] * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i,p,q], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β) * J[q]; limits = [0, β / 2], dims = dims)
        v[i,p,q], w[i,p,q], E[i,p,q] = variation((v, w) -> β == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper)
        v_guess = v_guesses == false ? v[i,p,q] : v_guesses
        w_guess = w_guesses == false ? w[i,p,q] : w_guesses
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J] | β = $(Ω[k]) [$k/$num_Ω]") end
            Σ[i,p,q,k] = holstein_memory(Ω[k], g[i,p,q], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[i,p,q], w[i,p,q]) * J[q] : polaron_propagator(t, v[i,p,q], w[i,p,q], β) * J[q]; dims = dims) * J[q]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::Number, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_ω, num_β, num_Ω = length(α), length(ω), length(β), length(Ω)

    if verbose N, n = num_α * num_ω * num_β * num_Ω, 1 end

    g = Matrix{Any}(undef, num_α, num_ω)
    v = Array{Any}(undef, num_α, num_ω, num_β)
    w = Array{Any}(undef, num_α, num_ω, num_β)
    E = Array{Any}(undef, num_α, num_ω, num_β)
    Σ = Array{Any}(undef, num_α, num_ω, num_β, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for i in eachindex(α), p in eachindex(ω)
        g[i,p] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J * E0_pu, dims))
        for j in eachindex(β)
            S(v, w) = holstein_S(v, w, g[i,p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J : polaron_propagator(τ, v, w, β[j]) * J; limits = [0, β[j] / 2], dims = dims)
            v[i,p,j], w[i,p,j], E[i,p,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            v_guess = v_guesses == false ? v[i,p,j] : v_guesses
            w_guess = w_guesses == false ? w[i,p,j] : w_guesses
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[i,p,j,k] = holstein_memory(Ω[k], g[i,p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,p,j], w[i,p,j]) * J : polaron_propagator(t, v[i,p,j], w[i,p,j], β[j]) * J; dims = dims) * J
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::Number, J::AbstractArray, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_J, num_β, num_Ω = length(α), length(J), length(β), length(Ω)

    if verbose N, n = num_α * num_J * num_β * num_Ω, 1 end

    g = Matrix{Any}(undef, num_α, num_J)
    v = Array{Any}(undef, num_α, num_J, num_β)
    w = Array{Any}(undef, num_α, num_J, num_β)
    E = Array{Any}(undef, num_α, num_J, num_β)
    Σ = Array{Any}(undef, num_α, num_J, num_β, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω + 1 / β[1] : w_guesses

    for i in eachindex(α), q in eachindex(J)
        g[i,q] = pustrip.(holstein_coupling(α[i], ω * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            S(v, w) = holstein_S(v, w, g[i,q], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[i,q,j], w[i,q,j], E[i,q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            v_guess = v_guesses == false ? v[i,q,j] : v_guesses
            w_guess = w_guesses == false ? w[i,q,j] : w_guesses
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | J = $(J[q]) [$q/$num_J] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[i,q,j,k] = holstein_memory(Ω[k], g[i,q], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,q,j], w[i,q,j]) * J[q] : polaron_propagator(t, v[i,q,j], w[i,q,j], β[j]) * J[q]; dims = dims) * J[q]
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Holstein(ω * ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::AbstractArray, J::AbstractArray, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_ω, num_J, num_β, num_Ω = length(ω), length(J), length(β), length(Ω)

    if verbose N, n = num_ω * num_J * num_β * num_Ω, 1 end

    g = Matrix{Any}(undef, num_ω, num_J)
    v = Array{Any}(undef, num_ω, num_J, num_β)
    w = Array{Any}(undef, num_ω, num_J, num_β)
    E = Array{Any}(undef, num_ω, num_J, num_β)
    Σ = Array{Any}(undef, num_ω, num_J, num_β, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for p in eachindex(ω), q in eachindex(J)
        g[p,q] = pustrip.(holstein_coupling(α, ω[p] * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            S(v, w) = holstein_S(v, w, g[p,q], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[p,q,j], w[p,q,j], E[p,q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            v_guess = v_guesses == false ? v[p,q,j] : v_guesses
            w_guess = w_guesses == false ? w[p,q,j] : w_guesses
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[p,q,j,k] = holstein_memory(Ω[k], g[p,q], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[p,q,j], w[p,q,j]) * J[q] : polaron_propagator(t, v[p,q,j], w[p,q,j], β[j]) * J[q]; dims = dims) * J[q]
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Holstein(ω .* ω0_pu, J * E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::AbstractArray, ω::AbstractArray, J::AbstractArray, β::AbstractArray, Ω::AbstractArray; dims = 3, v_guesses = false, w_guesses = false, upper = Inf, verbose = false)
    ω = pustrip.(ω)
    J = pustrip.(J)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_ω, num_J, num_β, num_Ω = length(α), length(ω), length(J), length(β), length(Ω)

    if verbose N, n = num_α * num_ω * num_J * num_β * num_Ω, 1 end

    g = Array{Any}(undef, num_α, num_ω, num_J)
    v = Array{Any}(undef, num_α, num_ω, num_J, num_β)
    w = Array{Any}(undef, num_α, num_ω, num_J, num_β)
    E = Array{Any}(undef, num_α, num_ω, num_J, num_β)
    Σ = Array{Any}(undef, num_α, num_ω, num_J, num_β, num_Ω)

    v_guess = v_guesses == false ? 2 * dims + 1 / β[1] : v_guesses
    w_guess = w_guesses == false ? ω[1] + 1 / β[1] : w_guesses

    for i in eachindex(α), p in eachindex(ω), q in eachindex(J)
        g[i,p,q] = pustrip(holstein_coupling(α[i], ω[p] * ω0_pu, J[q] * E0_pu, dims))
        for j in eachindex(β)
            S(v, w) = holstein_S(v, w, g[i,p,q], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * J[q] : polaron_propagator(τ, v, w, β[j]) * J[q]; limits = [0, β[j] / 2], dims = dims)
            v[i,p,q,j], w[i,p,q,j], E[i,p,q,j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J[q] + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J[q] + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper)
            v_guess = v_guesses == false ? v[i,p,q,j] : v_guesses
            w_guess = w_guesses == false ? w[i,p,q,j] : w_guesses
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | J = $(J[q]) [$q/$num_J] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[i,p,q,j,k] = holstein_memory(Ω[k], g[i,p,q], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,p,q,j], w[i,p,q,j]) * J[q] : polaron_propagator(t, v[i,p,q,j], w[i,p,q,j], β[j]) * J[q]; dims = dims) * J[q]
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
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
