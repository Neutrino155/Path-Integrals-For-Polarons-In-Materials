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

function holstein(α::Number, ω::Number, J::Number; dims = 3, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = holstein_S(v, w, g, τ -> phonon_propagator(τ, ω), (τ, v, w) -> polaron_propagator(τ, v, w); dims = dims)
    v, w, E = variation((v, w) -> - 2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3), 2 * dims + ω, ω)
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, Inf / E0_pu, 0 * ω0_pu, 0 * ω0_pu)
end

function holstein(α, ω::Number, J::Number; dims = 3, verbose = false)
    num_α = length(α)
    ω = pustrip(ω)
    J = pustrip(J)
    g = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    for i in 1:num_α
        g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> phonon_propagator(τ, ω), (τ, v, w) -> polaron_propagator(τ, v, w); dims = dims)
        v[i], w[i], E[i] = variation((v, w) -> -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3), i > 1 ? v[i-1] : 2 * dims, i > 1 ? w[i-1] : ω)
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, Inf / E0_pu, 0 * ω0_pu, 0 * ω0_pu)
end

function holstein(α::Number, ω, J::Number; dims = 3, verbose = false)
    num_ω = length(ω)
    ω = pustrip.(ω)
    J = pustrip(J)
    g = Vector{Any}(undef, num_ω)
    v = Vector{Any}(undef, num_ω)
    w = Vector{Any}(undef, num_ω)
    E = Vector{Any}(undef, num_ω)
    for i in 1:num_ω
        g[i] = pustrip(holstein_coupling(α, ω[i] * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> phonon_propagator(τ, ω[i]), (τ, v, w) -> polaron_propagator(τ, v, w); dims = dims)
        v[i], w[i], E[i] = variation((v, w) -> -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3), i > 1 ? v[i-1] : 2 * dims, i > 1 ? w[i-1] : ω[i])
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, Inf / E0_pu, 0 * ω0_pu, 0 * ω0_pu)
end

function holstein(α, ω, J::Number; dims = 3, verbose = false)
    num_α = length(α)
    num_ω = length(ω)
    ω = pustrip.(ω)
    J = pustrip(J)
    g = Matrix{Any}(undef, num_α, num_ω)
    v = Matrix{Any}(undef, num_α, num_ω)
    w = Matrix{Any}(undef, num_α, num_ω)
    E = Matrix{Any}(undef, num_α, num_ω)
    for i in 1:num_α, j in 1:num_ω
        g[i, j] = pustrip(holstein_coupling(α[i], ω[j] * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i, j], τ -> phonon_propagator(τ, ω[j]), (τ, v, w) -> polaron_propagator(τ, v, w), dims = dims)
        v[i, j], w[i, j], E[i, j] = variation((v, w) -> -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3), i > 1 ? v[i-1, j] : 2 * dims, i > 1 ? w[i-1, j] : ω[j])
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, Inf / E0_pu, 0 * ω0_pu, 0 * ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β::Number; dims = 3, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = holstein_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β); limits = [0, β / 2], dims = dims)
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), 2 * dims + 1 / β, ω + 1 / β)
    Σ = holstein_memory(eps(), g, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) : polaron_propagator(t, v, w, β), dims = dims) * ω
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function holstein(α, ω::Number, J::Number, β::Number; dims = 3, verbose = false)
    num_α = length(α)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    g = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Vector{Any}(undef, num_α)
    @inbounds for i in 1:num_α
        @fastmath g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β), limits = [0, β/2]; dims = dims)
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), i > 1 ? v[i-1] : 2 * dims + 1 / β, i > 1 ? w[i-1] : ω + 1 / β)
        @fastmath Σ[i] = holstein_memory(eps(), g[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) : polaron_propagator(t, v[i], w[i], β); dims = dims) * ω
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β; dims = 3, verbose = false)
    num_β = length(β)
    if verbose N, n = num_β, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)
    @inbounds for i in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[i]) [$i/$num_β]") end
        S(v, w) = holstein_S(v, w, g, τ -> β[i] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[i]), (τ, v, w) -> β[i] == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β[i]), limits = [0, β[i]/2]; dims = dims)
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β[i] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[i]) * dims / 3 - (S(v, w) - S₀(v, w, β[i]) * dims / 3), i > 1 ? v[i-1] : 2 * dims + 1 / β[i], i > 1 ? w[i-1] : ω + 1 / β[i])
        @fastmath Σ[i] = holstein_memory(eps(), g, t -> β[i] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[i]), t -> β[i] == Inf ? polaron_propagator(t, v[i], w[i]) : polaron_propagator(t, v[i], w[i], β[i]); dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function holstein(α, ω::Number, J::Number, β; dims = 3, verbose = false)
    num_α = length(α)
    num_β = length(β)
    if verbose N, n = num_α * num_β, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    g = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Matrix{Any}(undef, num_α, num_β)
    @inbounds for j in 1:num_β, i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
        @fastmath g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β[j]), limits = [0, β[j]/2]; dims = dims)
        @fastmath v[i, j], w[i, j], E[i, j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), i > 1 ? v[i-1, j] : 2 * dims + 1 / β[j], i > 1 ? w[i-1, j] : ω + 1 / β[j])
        @fastmath Σ[i, j] = holstein_memory(eps(), g[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i, j], w[i, j]) : polaron_propagator(t, v[i, j], w[i, j], β[j]); dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β::Number, Ω::Number; dims = 3, verbose = false)
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip(β * ħ_pu / E0_pu) 
    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))
    S(v, w) = holstein_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β); limits = [0, β / 2], dims = dims)
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), 2 * dims, ω)
    Σ = holstein_memory(Ω, g, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) : polaron_propagator(t, v, w, β); dims = dims) * ω
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α, ω::Number, J::Number, β::Number, Ω::Number; dims = 3, verbose = false)
    num_α = length(α)
    if verbose N, n = num_α, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    g = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Vector{Any}(undef, num_α)
    @inbounds for i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α]") end
        @fastmath g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β), limits = [0, β/2]; dims = dims)
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), i > 1 ? v[i-1] : 2 * dims + 1 / β, i > 1 ? w[i-1] : ω + 1 / β)
        @fastmath Σ[i] = holstein_memory(Ω, g[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) : polaron_propagator(t, v[i], w[i], β); dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β, Ω::Number; dims = 3, verbose = false)
    num_β = length(β)
    if verbose N, n = num_β, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)

    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))

    @inbounds for j in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = holstein_S(v, w, g, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β[j]), limits = [0, β[j]/2]; dims = dims)
        @fastmath v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), i > 1 ? v[j-1] : 2 * dims + 1 / β[j], i > 1 ? w[j-1] : ω + 1 / β[j])
        @fastmath Σ[j] = holstein_memory(Ω, g, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) : polaron_propagator(t, v[j], w[j], β[j]); dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β::Number, Ω; dims = 3, verbose = false)
    num_Ω = length(Ω)
    if verbose N, n = num_Ω, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Σ = Vector{Any}(undef, num_Ω)

    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))

    S(v, w) = holstein_S(v, w, g, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β), limits = [0, β/2]; dims = dims)
    v, w, E = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), 2 * dims + 1 / β, ω + 1 / β)
    @inbounds for k in 1:num_Ω
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | Ω = $(Ω[k]) [$k/$num_Ω]") end
        @fastmath Σ[k] = holstein_memory(Ω[k], g, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) : polaron_propagator(t, v, w, β); dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α, ω::Number, J::Number, β::Number, Ω; dims = 3, verbose = false)
    num_α = length(α)
    num_Ω = length(Ω)
    if verbose N, n = num_α * num_Ω, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    g = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Matrix{Any}(undef, num_α, num_Ω)
    @inbounds for i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α]") end
        @fastmath g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β), limits = [0, β/2]; dims = dims)
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), i > 1 ? v[i-1] : 2 * dims + 1 / β, i > 1 ? w[i-1] : ω + 1 / β)
        if verbose print("\e[1F") end
        @inbounds for k in 1:num_Ω
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            @fastmath Σ[i, k] = holstein_memory(Ω[k], g[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) : polaron_propagator(t, v[i], w[i], β); dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end

    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α::Number, ω::Number, J::Number, β, Ω; dims = 3, verbose = false)
    num_β = length(β)
    num_Ω = length(Ω)
    if verbose N, n = num_β * num_Ω, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Matrix{Any}(undef, num_β, num_Ω)

    g = pustrip(holstein_coupling(α, ω * ω0_pu, J * E0_pu, dims))

    @inbounds for j in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = holstein_S(v, w, g, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β[j]), limits = [0, β[j]/2]; dims = dims)
        @fastmath v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), j > 1 ? v[j-1] : 2 * dims + 1 / β[j], j > 1 ? w[j-1] : ω + 1 / β[j])
        if verbose print("\e[1F") end
        @inbounds for k in 1:num_Ω
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            @fastmath Σ[j, k] = holstein_memory(Ω[k], g, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) : polaron_propagator(t, v[j], w[j], β[j]); dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α, ω::Number, J::Number, β, Ω::Number; dims = 3, verbose = false)
    num_α = length(α)
    num_β = length(β)
    if verbose N, n = num_α * num_β, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    g = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Matrix{Any}(undef, num_α, num_β)
    @inbounds for j in 1:num_β, i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
        @fastmath g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β[j]), limits = [0, β[j]/2]; dims = dims)
        @fastmath v[i, j], w[i, j], E[i, j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), i > 1 && j > 1 ? v[i-1, j-1] : 2 * dims + 1 / β[j], i > 1 && j > 1 ? w[i-1, j-1] : ω + 1 / β[j])
        @fastmath Σ[i, j] = holstein_memory(Ω, g[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i, j], w[i, j]) : polaron_propagator(t, v[i, j], w[i, j], β[j]); dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(α, ω::Number, J::Number, β, Ω; dims = 3, verbose = false)
    num_α = length(α)
    num_β = length(β)
    num_Ω = length(Ω)
    if verbose N, n = num_α * num_β * num_Ω, 1 end
    ω = pustrip(ω)
    J = pustrip(J)
    β = pustrip.(β .* ħ_pu / E0_pu)
    g = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Array{Any}(undef, num_α, num_β, num_Ω)
    @inbounds for j in 1:num_β, i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
        @fastmath g[i] = pustrip(holstein_coupling(α[i], ω * ω0_pu, J * E0_pu, dims))
        S(v, w) = holstein_S(v, w, g[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β[j]), limits = [0, β[j]/2]; dims = dims)
        @fastmath v[i, j], w[i, j], E[i, j] = variation((v, w) -> β[j] == Inf ? -2 * dims * J + E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : -2 * dims * J + E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), i > 1 && j > 1 ? v[i-1, j-1] : 2 * dims + 1 / β[j], i > 1 && j > 1 ? w[i-1, j-1] : ω + 1 / β[j])
        if verbose print("\e[1F") end
        @inbounds for k in 1:num_Ω
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            @fastmath Σ[i, j, k] = holstein_memory(Ω[k], g[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i, j], w[i, j]) : polaron_propagator(t, v[i, j], w[i, j], β[j]); dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Holstein(ω .* ω0_pu, J .* E0_pu, dims, g .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function holstein(material::Material; verbose = false)
    α = holstein_alpha.(material.g, material.ω_LO, material.J, material.d)
    return holstein(α, material.ω_LO, material.J; dims = material.d, verbose = verbose)
end

function holstein(material::Material, T; verbose = false)
    α = holstein_alpha.(material.g, material.ω_LO, material.J, material.d)
    return holstein(α, material.ω_LO, material.J, pustrip.(1 ./ (kB_pu .* T)); dims = material.d, verbose = verbose)
end

function holstein(material::Material, T, Ω; verbose = false)
    α = holstein_alpha.(material.g, material.ω_LO, material.J, material.d)
    return holstein(α, material.ω_LO, material.J, pustrip.(1 ./ (kB_pu .* T)), pustrip.(Ω); dims = material.d, verbose = verbose)
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
    return norm(coupling)^2 / (2π)^dims * integral
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
