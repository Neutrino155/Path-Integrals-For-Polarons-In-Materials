# Frohlich.jl

mutable struct Frohlich

    # Hamiltonian parameters
    ω
    Mₖ
    α

    # Polaron properties
    E
    v
    w
    β
    Ω
    Σ

    function Frohlich(x...)
        new(reduce_array.(x)...)
    end
end

frohlich(α::Number; dims = 3) = frohlich(α, 1; dims = dims, verbose = false)

frohlich(α::AbstractArray; dims = 3) = frohlich(α, 1; dims = dims, verbose = false)

function frohlich(α::Number, ω::Number; mb = 1m0_pu, dims = 3, verbose = false)
    ω = pustrip(ω)
    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = frohlich_S(v, w, Mₖ, τ -> phonon_propagator(τ, ω), (τ, v, w) -> polaron_propagator(τ, v, w) * ω; dims = dims)
    v, w, E = variation((v, w) -> E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3), α < 7 ? 3 + α / 4 : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4, 2 + tanh((6 - α) / 3))
    return Frohlich(ω * ω0_pu, Mₖ * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, Inf / E0_pu, 0 * ω0_pu, 0 * ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number; mb = 1m0_pu, dims = 3, verbose = false)
    num_α = length(α)
    ω = pustrip(ω)
    Mₖ = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    for i in 1:num_α
        Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> phonon_propagator(τ, ω), (τ, v, w) -> polaron_propagator(τ, v, w) * ω; dims = dims)
        v[i], w[i], E[i] = variation((v, w) -> E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3), α[i] < 7 ? 3 + α[i] / 4 : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4, 2 + tanh((6 - α[i]) / 3))
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, Inf / E0_pu, 0 * ω0_pu, 0 * ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    @assert length(α) == length(ω)
    ω = pustrip.(ω)
    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))
    S(v, w) = sum(frohlich_S(v, w, Mₖ[j], τ -> phonon_propagator(τ, ω[j]), (τ, v, w) -> polaron_propagator(τ, v, w) * ω[j]; dims = dims) for j in eachindex(ω))
    v, w, E = variation((v, w) -> E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3), sum(α) < 7 ? 3 + sum(α) / 4 : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4, 2 + tanh((6 - sum(α)) / 3))
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, Inf / E0_pu, 0 * ω0_pu, 0 * ω0_pu)
end

function frohlich(α::Number, ω::Number, β::Number; mb = 1m0_pu, dims = 3, verbose = false)
    ω = pustrip(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = frohlich_S(v, w, Mₖ, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, β / 2])
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - α) / 3) + 1 / β)
    Σ = frohlich_memory(eps(), Mₖ, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * ω : polaron_propagator(t, v, w, β) * ω; dims = dims) * ω
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::Number; mb = 1m0_pu, dims = 3, verbose = false)
    @assert length(α) == length(ω)
    ω = pustrip.(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))
    S(v, w) = sum(frohlich_S(v, w, Mₖ[j], τ -> β == Inf ? phonon_propagator(τ, ω[j]) : phonon_propagator(τ, ω[j], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω[j] : polaron_propagator(τ, v, w, β) * ω[j]; dims = dims, limits = [0, β / 2]) for j in eachindex(ω))
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - sum(α)) / 3) + 1 / β)
    Σ = sum(frohlich_memory(eps(), Mₖ[j], t -> β == Inf ? phonon_propagator(t, ω[j]) : phonon_propagator(t, ω[j], β), t -> β == Inf ? polaron_propagator(t, v, w) * ω[j] : polaron_propagator(t, v, w, β) * ω[j]) * ω[j] for j in eachindex(ω); dims = dims)
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::Number; mb = 1m0_pu, dims = 3, verbose = false)
    num_α = length(α)
    ω = pustrip(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Vector{Any}(undef, num_α)
    @inbounds for i in 1:num_α
        @fastmath Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, β/2])
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), α[i] < 7 ? 3 + α[i] / 4 + 1 / β : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - α[i]) / 3) + 1 / β)
        @fastmath Σ[i] = frohlich_memory(eps(), Mₖ[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) * ω : polaron_propagator(t, v[i], w[i], β) * ω; dims = dims) * ω
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    num_β = length(β)
    if verbose N, n = num_β, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)
    @inbounds for i in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[i]) [$i/$num_β]") end
        S(v, w) = frohlich_S(v, w, Mₖ, τ -> β[i] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[i]), (τ, v, w) -> β[i] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[i]) * ω, dims = dims, limits = [0, β[i]/2])
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β[i] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[i]) - (S(v, w) - S₀(v, w, β[i])), α < 7 ? 3 + α / 4 + 1 / β[i] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[i], 2 + tanh((6 - α) / 3) + 1 / β[i])
        @fastmath Σ[i] = frohlich_memory(eps(), Mₖ, t -> β[i] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[i]), t -> β[i] == Inf ? polaron_propagator(t, v[i], w[i]) * ω : polaron_propagator(t, v[i], w[i], β[i]) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ * E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    @assert length(α) == length(ω)
    num_β = length(β)
    if verbose N, n = num_β, 1 end
    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)
    @inbounds for i in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[i]) [$i/$num_β]") end
        S(v, w) = sum(frohlich_S(v, w, Mₖ[j], τ -> β[i] == Inf ? phonon_propagator(τ, ω[j]) : phonon_propagator(τ, ω[j], β[i]), (τ, v, w) -> β[i] == Inf ? polaron_propagator(τ, v, w) * ω[j] : polaron_propagator(τ, v, w, β[i]) * ω[j]; dims = dims, limits = [0, β[i]/2]) for j in eachindex(ω))
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β[i] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[i]) - (S(v, w) - S₀(v, w, β[i])), sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β[i] : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[i], 2 + tanh((6 - sum(α)) / 3) + 1 / β[i])
        @fastmath Σ[i] = sum(frohlich_memory(eps(), Mₖ[j], t -> β[i] == Inf ? phonon_propagator(t, ω[j]) : phonon_propagator(t, ω[j], β[i]), t -> β[i] == Inf ? polaron_propagator(t, v[i], w[i]) * ω[j] : polaron_propagator(t, v[i], w[i], β[i]) * ω[j]; dims = dims) * ω[j] for j in eachindex(ω))
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    num_α = length(α)
    num_β = length(β)
    if verbose N, n = num_α * num_β, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Mₖ = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Matrix{Any}(undef, num_α, num_β)
    @inbounds for j in 1:num_β, i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
        @fastmath Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, β[j]/2])
        @fastmath v[i, j], w[i, j], E[i, j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) - (S(v, w) - S₀(v, w, β[j])), α[i] < 7 ? 3 + α[i] / 4 + 1 / β[j] : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j], 2 + tanh((6 - α[i]) / 3) + 1 / β[j])
        @fastmath Σ[i, j] = frohlich_memory(eps(), Mₖ[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i, j], w[i, j]) * ω : polaron_propagator(t, v[i, j], w[i, j], β[j]) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, 0 * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::Number, Ω::Number, mb = 1m0_pu, dims = 3, verbose = false)
    ω = pustrip(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = frohlich_S(v, w, Mₖ, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, β / 2])
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - α) / 3) + 1 / β)
    Σ = frohlich_memory(Ω, Mₖ, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * ω : polaron_propagator(t, v, w, β) * ω; dims = dims) * ω
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, verbose = false)
    @assert length(α) == length(ω)
    ω = pustrip.(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))
    S(v, w) = sum(frohlich_S.(v, w, Mₖ, τ -> β == Inf ? phonon_propagator.(τ, ω) : phonon_propagator.(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) .* ω : polaron_propagator(τ, v, w, β) .* ω; dims = dims, limits = [0, β / 2]))
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - sum(α)) / 3) + 1 / β)
    Σ = sum(frohlich_memory.(Ω, Mₖ, t -> β == Inf ? phonon_propagator.(t, ω) : phonon_propagator.(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) .* ω : polaron_propagator(t, v, w, β) .* ω; dims = dims) .* ω)
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, verbose = false)
    num_α = length(α)
    if verbose N, n = num_α, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Mₖ = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Vector{Any}(undef, num_α)
    @inbounds for i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α]") end
        @fastmath Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, β/2])
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), α[i] < 7 ? 3 + α[i] / 4 + 1 / β : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - α[i]) / 3) + 1 / β)
        @fastmath Σ[i] = frohlich_memory(Ω, Mₖ[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) * ω : polaron_propagator(t, v[i], w[i], β) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::AbstractArray, Ω::Number; mb = 1m0_pu, dims = 3, verbose = false)
    num_β = length(β)
    if verbose N, n = num_β, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)

    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))

    @inbounds for j in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = frohlich_S(v, w, Mₖ, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, β[j]/2])
        @fastmath v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) - (S(v, w) - S₀(v, w, β[j])), α < 7 ? 3 + α / 4 + 1 / β[j] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j], 2 + tanh((6 - α) / 3) + 1 / β[j])
        @fastmath Σ[j] = frohlich_memory(Ω, Mₖ, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) * ω : polaron_propagator(t, v[j], w[j], β[j]) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::AbstractArray, Ω::Number; mb = 1m0_pu, dims = 3, verbose = false)
    @assert length(α) == length(ω)
    num_β = length(β)
    if verbose N, n = num_β, 1 end
    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)

    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))

    @inbounds for j in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = sum(frohlich_S.(v, w, Mₖ, τ -> β[j] == Inf ? phonon_propagator.(τ, ω) : phonon_propagator.(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) .* ω : polaron_propagator(τ, v, w, β[j]) .* ω; dims = dims, limits = [0, β[j]/2]))
        @fastmath v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) - (S(v, w) - S₀(v, w, β[j])), sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β[j] : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j], 2 + tanh((6 - sum(α)) / 3) + 1 / β[j])
        @fastmath Σ[j] = sum(frohlich_memory.(Ω, Mₖ, t -> β[j] == Inf ? phonon_propagator.(t, ω) : phonon_propagator.(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) .* ω : polaron_propagator(t, v[j], w[j], β[j]) .* ω; dims = dims) .* ω)
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::Number, Ω::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    num_Ω = length(Ω)
    if verbose N, n = num_Ω, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Σ = Vector{Any}(undef, num_Ω)

    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))

    S(v, w) = frohlich_S(v, w, Mₖ, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, β/2])
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - α) / 3) + 1 / β)
    @inbounds for k in 1:num_Ω
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | Ω = $(Ω[k]) [$k/$num_Ω]") end
        @fastmath Σ[k] = frohlich_memory(Ω[k], Mₖ, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * ω : polaron_propagator(t, v, w, β) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::Number, Ω::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    @assert length(α) == length(ω)
    num_Ω = length(Ω)
    if verbose N, n = num_Ω, 1 end
    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Σ = Vector{Any}(undef, num_Ω)

    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))

    S(v, w) = sum(frohlich_S.(v, w, Mₖ, τ -> β == Inf ? phonon_propagator.(τ, ω) : phonon_propagator.(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) .* ω : polaron_propagator(τ, v, w, β) .* ω; dims = dims, limits = [0, β/2]))
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - sum(α)) / 3) + 1 / β)
    @inbounds for k in 1:num_Ω
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | Ω = $(Ω[k]) [$k/$num_Ω]") end
        @fastmath Σ[k] = sum(frohlich_memory.(Ω[k], Mₖ, t -> β == Inf ? phonon_propagator.(t, ω) : phonon_propagator.(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) .* ω : polaron_propagator(t, v, w, β) .* ω; dims = dims) .* ω)
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::Number, Ω::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    num_α = length(α)
    num_Ω = length(Ω)
    if verbose N, n = num_α * num_Ω, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Mₖ = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Matrix{Any}(undef, num_α, num_Ω)
    @inbounds for i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α]") end
        @fastmath Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, β/2])
        @fastmath v[i], w[i], E[i] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), α[i] < 7 ? 3 + α[i] / 4 + 1 / β : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β, 2 + tanh((6 - α[i]) / 3) + 1 / β)
        if verbose print("\e[1F") end
        @inbounds for k in 1:num_Ω
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            @fastmath Σ[i, k] = frohlich_memory(Ω[k], Mₖ[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) * ω : polaron_propagator(t, v[i], w[i], β) * ω; dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end

    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::AbstractArray, Ω::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    num_β = length(β)
    num_Ω = length(Ω)
    if verbose N, n = num_β * num_Ω, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Matrix{Any}(undef, num_β, num_Ω)

    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))

    @inbounds for j in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = frohlich_S(v, w, Mₖ, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, β[j]/2])
        @fastmath v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) - (S(v, w) - S₀(v, w, β[j])), α < 7 ? 3 + α / 4 + 1 / β[j] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j], 2 + tanh((6 - α) / 3) + 1 / β[j])
        if verbose print("\e[1F") end
        @inbounds for k in 1:num_Ω
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            @fastmath Σ[j, k] = frohlich_memory(Ω[k], Mₖ, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) * ω : polaron_propagator(t, v[j], w[j], β[j]) * ω; dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::AbstractArray, Ω::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    @assert length(α) == length(ω)
    num_β = length(β)
    num_Ω = length(Ω)
    if verbose N, n = num_β * num_Ω, 1 end
    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Matrix{Any}(undef, num_β, num_Ω)

    Mₖ = pustrip.(frohlich_coupling.(α, ω .* ω0_pu, mb; dims = dims))

    @inbounds for j in 1:num_β
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = sum(frohlich_S.(v, w, Mₖ, τ -> β[j] == Inf ? phonon_propagator.(τ, ω) : phonon_propagator.(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) .* ω : polaron_propagator(τ, v, w, β[j]) .* ω; dims = dims, limits = [0, β[j]/2]))
        @fastmath v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) - (S(v, w) - S₀(v, w, β[j])), sum(α) < 7 ? 3 + sum(α) / 4 + 1 / β[j] : 4 * sum(α)^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j], 2 + tanh((6 - sum(α)) / 3) + 1 / β[j])
        if verbose print("\e[1F") end
        @inbounds for k in 1:num_Ω
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            @fastmath Σ[j, k] = sum(frohlich_memory.(Ω[k], Mₖ, t -> β[j] == Inf ? phonon_propagator.(t, ω) : phonon_propagator.(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) .* ω : polaron_propagator(t, v[j], w[j], β[j]) .* ω; dims = dims) .* ω)
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::AbstractArray, Ω::Number; mb = 1m0_pu, dims = 3, verbose = false)
    num_α = length(α)
    num_β = length(β)
    if verbose N, n = num_α * num_β, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Mₖ = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Matrix{Any}(undef, num_α, num_β)
    @inbounds for j in 1:num_β, i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
        @fastmath Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, β[j]/2])
        @fastmath v[i, j], w[i, j], E[i, j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) - (S(v, w) - S₀(v, w, β[j])), α[i] < 7 ? 3 + α[i] / 4 + 1 / β[j] : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j], 2 + tanh((6 - α[i]) / 3) + 1 / β[j])
        @fastmath Σ[i, j] = frohlich_memory(Ω, Mₖ[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i, j], w[i, j]) * ω : polaron_propagator(t, v[i, j], w[i, j], β[j]) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::AbstractArray, Ω::AbstractArray; mb = 1m0_pu, dims = 3, verbose = false)
    num_α = length(α)
    num_β = length(β)
    num_Ω = length(Ω)
    if verbose N, n = num_α * num_β * num_Ω, 1 end
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu)
    Mₖ = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Array{Any}(undef, num_α, num_β, num_Ω)
    @inbounds for j in 1:num_β, i in 1:num_α
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
        @fastmath Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, β[j]/2])
        @fastmath v[i, j], w[i, j], E[i, j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) - (S(v, w) - S₀(v, w, β[j])), α[i] < 7 ? 3 + α[i] / 4 + 1 / β[j] : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j], 2 + tanh((6 - α[i]) / 3) + 1 / β[j])
        if verbose print("\e[1F") end
        @inbounds for k in 1:num_Ω
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            @fastmath Σ[i, j, k] = frohlich_memory(Ω[k], Mₖ[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i, j], w[i, j]) * ω : polaron_propagator(t, v[i, j], w[i, j], β[j]) * ω; dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(material::Material; dims = 3, verbose = false)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO; dims = dims, verbose = verbose)
end

function frohlich(material::Material, T; dims = 3, verbose = false)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO, pustrip.(1 ./ (kB_pu .* T)); dims = dims, verbose = verbose)
end

function frohlich(material::Material, T, Ω; dims = 3, verbose = false)
    α = frohlich_alpha.(material.ϵ_optic, material.ϵ_total, material.ϵ_ionic, material.ω_LO, material.mb)
    return frohlich(α, material.ω_LO, pustrip.(1 ./ (kB_pu .* T)), pustrip.(Ω); dims = dims, verbose = verbose)
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
