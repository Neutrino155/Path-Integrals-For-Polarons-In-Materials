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

function frohlich(α::Number, ω::Number, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = frohlich_S(v, w, Mₖ, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, sqrt(β / 2)])
    v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
    w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
    Σ = isinf(β) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * ω : polaron_propagator(t, v, w, β) * ω; dims = dims) * ω

    return Frohlich(ω * ω0_pu, Mₖ * E0_pu, α, E * E0_pu, v * ω0_pu, w * ω0_pu, β / E0_pu, Ω * ω0_pu, Σ * ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α = length(α)
    if verbose N, n = num_α, 1 end
    Mₖ = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Vector{Any}(undef, num_α)

    for i in eachindex(α)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α]") end
        Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, sqrt(β / 2)])
        v_guess = v_guesses == false ? α[i] < 7 ? 3 + α[i] / 4 + 1 / β : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α[i]) / 3) + 1 / β : w_guesses
        v[i], w[i], E[i] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        Σ[i] = isinf(β) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) * ω : polaron_propagator(t, v[i], w[i], β) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::AbstractArray, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip.(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_ω = length(ω)
    if verbose N, n = num_ω, 1 end

    Mₖ = Vector{Any}(undef, num_ω)
    v = Vector{Any}(undef, num_ω)
    w = Vector{Any}(undef, num_ω)
    E = Vector{Any}(undef, num_ω)
    Σ = Vector{Any}(undef, num_ω)

    for p in eachindex(ω)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω]") end
        Mₖ[p] = pustrip(frohlich_coupling(α, ω[p] * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β) * ω[p]; dims = dims, limits = [0, sqrt(β / 2)])
        v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β : w_guesses
        v[p], w[p], E[p] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        Σ[p] = isinf(β) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ[p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[p], w[p]) * ω[p] : polaron_propagator(t, v[p], w[p], β) * ω[p]; dims = dims) * ω[p]
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::AbstractArray, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_β = length(β)
    if verbose N, n = num_β, 1 end

    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Vector{Any}(undef, num_β)

    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    for j in eachindex(β)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β]") end
        S(v, w) = frohlich_S(v, w, Mₖ, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, sqrt(β[j] / 2)])
        v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β[j] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β[j] : w_guesses
        v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        Σ[j] = isinf(β[j]) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) * ω : polaron_propagator(t, v[j], w[j], β[j]) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ * E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::Number, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_Ω = length(Ω)
    if verbose N, n = num_Ω, 1 end

    Σ = Vector{Any}(undef, num_Ω)

    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    S(v, w) = frohlich_S(v, w, Mₖ, τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, sqrt(β / 2)])
    v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
    w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β : w_guesses
    v, w, E = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
    for k in eachindex(Ω)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | Ω = $(Ω[k]) [$k/$num_Ω]") end
        Σ[k] = isinf(β) && Ω[k] == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ, t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v, w) * ω : polaron_propagator(t, v, w, β) * ω; dims = dims) * ω
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::Number, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip.(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_ω = length(α), length(ω)
    if verbose N, n = num_α * num_ω, 1 end

    Mₖ = Matrix{Any}(undef, num_α, num_ω)
    v = Matrix{Any}(undef, num_α, num_ω)
    w = Matrix{Any}(undef, num_α, num_ω)
    E = Matrix{Any}(undef, num_α, num_ω)
    Σ = Matrix{Any}(undef, num_α, num_ω)

    for i in eachindex(α), p in eachindex(ω)
        if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω]") end
        Mₖ[i,p] = pustrip(frohlich_coupling(α[i], ω[p] * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i,p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β) * ω[p]; dims = dims, limits = [0, sqrt(β / 2)])
        v_guess = v_guesses == false ? α[i] < 7 ? 3 + α[i] / 4 + 1 / β : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α[i]) / 3) + 1 / β : w_guesses
        v[i,p], w[i,p], E[i,p] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        Σ[i,p] = isinf(β) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ[i,p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[i,p], w[i,p]) * ω[p] : polaron_propagator(t, v[i,p], w[i,p], β) * ω[p]; dims = dims) * ω[p]
        if verbose n += 1; print("\e[1F") end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::AbstractArray, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_β = length(α), length(β)
    if verbose N, n = num_α * num_β, 1 end

    Mₖ = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Matrix{Any}(undef, num_α, num_β)

    for i in eachindex(α)
        Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, sqrt(β[j] / 2)])
            v_guess = v_guesses == false ? α[i] < 7 ? 3 + α[i] / 4 + 1 / β[j] : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
            w_guess = w_guesses == false ? 2 + tanh((6 - α[i]) / 3) + 1 / β[j] : w_guesses
            v[i,j], w[i,j], E[i,j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
            Σ[i,j] = isinf(β[j]) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,j], w[i,j]) * ω : polaron_propagator(t, v[i,j], w[i,j], β[j]) * ω; dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::Number, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_Ω = length(α), length(Ω)
    if verbose N, n = num_α * num_Ω, 1 end

    Mₖ = Vector{Any}(undef, num_α)
    v = Vector{Any}(undef, num_α)
    w = Vector{Any}(undef, num_α)
    E = Vector{Any}(undef, num_α)
    Σ = Matrix{Any}(undef, num_α, num_Ω)

    for i in eachindex(α)
        Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β) * ω; dims = dims, limits = [0, sqrt(β / 2)])
        v_guess = v_guesses == false ? α[i] < 7 ? 3 + α[i] / 4 + 1 / β : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α[i]) / 3) + 1 / β : w_guesses
        v[i], w[i], E[i] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[i,k] = isinf(β) && Ω[k] == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ[i], t -> β == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β), t -> β == Inf ? polaron_propagator(t, v[i], w[i]) * ω : polaron_propagator(t, v[i], w[i], β) * ω; dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::AbstractArray, β::AbstractArray, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_ω, num_β = length(ω), length(β)
    if verbose N, n = num_ω * num_β, 1 end

    Mₖ = Vector{Any}(undef, num_ω)
    v = Matrix{Any}(undef, num_ω, num_β)
    w = Matrix{Any}(undef, num_ω, num_β)
    E = Matrix{Any}(undef, num_ω, num_β)
    Σ = Matrix{Any}(undef, num_ω, num_β)

    for p in eachindex(ω)
        Mₖ[p] = pustrip(frohlich_coupling(α, ω[p] * ω0_pu, mb; dims = dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = frohlich_S(v, w, Mₖ[p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β[j]) * ω[p]; dims = dims, limits = [0, sqrt(β[j] / 2)])
            v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β[j] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
            w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β[j] : w_guesses
            v[p,j], w[p,j], E[p,j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
            Σ[p,j] = isinf(β[j]) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ[p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[p,j], w[p,j]) * ω[p] : polaron_propagator(t, v[p,j], w[p,j], β[j]) * ω[p]; dims = dims) * ω[p]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::AbstractArray, β::Number, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip.(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_ω, num_Ω = length(ω), length(Ω)
    if verbose N, n = num_ω * num_Ω, 1 end

    Mₖ = Vector{Any}(undef, num_ω)
    v = Vector{Any}(undef, num_ω)
    w = Vector{Any}(undef, num_ω)
    E = Vector{Any}(undef, num_ω)
    Σ = Matrix{Any}(undef, num_ω, num_Ω)

    for p in eachindex(ω)
        Mₖ[p] = pustrip(frohlich_coupling(α, ω[p] * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β) * ω[p]; dims = dims, limits = [0, sqrt(β / 2)])
        v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β : w_guesses
        v[p], w[p], E[p] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[p,k] = isinf(β) && Ω[k] == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ[p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[p], w[p]) * ω[p] : polaron_propagator(t, v[p], w[p], β) * ω[p]; dims = dims) * ω[p]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::Number, β::AbstractArray, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)
    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_β, num_Ω = length(β), length(Ω)
    if verbose N, n = num_β * num_Ω, 1 end

    v = Vector{Any}(undef, num_β)
    w = Vector{Any}(undef, num_β)
    E = Vector{Any}(undef, num_β)
    Σ = Matrix{Any}(undef, num_β, num_Ω)

    Mₖ = pustrip(frohlich_coupling(α, ω * ω0_pu, mb; dims = dims))
    for j in eachindex(β)
        S(v, w) = frohlich_S(v, w, Mₖ, τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, sqrt(β[j] / 2)])
        v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β[j] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β[j] : w_guesses
        v[j], w[j], E[j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[j,k] = isinf(β[j]) && Ω[k] == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ, t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[j], w[j]) * ω : polaron_propagator(t, v[j], w[j], β[j]) * ω; dims = dims) * ω
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω * ω0_pu, Mₖ * E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::Number, ω::AbstractArray, β::AbstractArray, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)

    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_ω, num_β, num_Ω = length(ω), length(β), length(Ω)
    if verbose N, n = num_ω * num_β * num_Ω, 1 end

    Mₖ = Vector{Any}(undef, num_ω)
    v = Matrix{Any}(undef, num_ω, num_β)
    w = Matrix{Any}(undef, num_ω, num_β)
    E = Matrix{Any}(undef, num_ω, num_β)
    Σ = Array{Any}(undef, num_ω, num_β, num_Ω)

    for p in eachindex(ω)
        Mₖ[p] = pustrip(frohlich_coupling(α, ω[p] * ω0_pu, mb; dims = dims))
        for j in eachindex(β)
            S(v, w) = frohlich_S(v, w, Mₖ[p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β[j]) * ω[p]; dims = dims, limits = [0, sqrt(β[j] / 2)])
            v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β[j] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
            w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β[j] : w_guesses
            v[p,j], w[p,j], E[p,j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[p,j,k] = isinf(β[j]) && Ω[k] == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ[p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[p,j], w[p,j]) * ω[p] : polaron_propagator(t, v[p,j], w[p,j], β[j]) * ω[p]; dims = dims) * ω[p]
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::Number, β::AbstractArray, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)

    ω = pustrip(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_β, num_Ω = length(α), length(β), length(Ω)
    if verbose N, n = num_α * num_β * num_Ω, 1 end

    Mₖ = Vector{Any}(undef, num_α)
    v = Matrix{Any}(undef, num_α, num_β)
    w = Matrix{Any}(undef, num_α, num_β)
    E = Matrix{Any}(undef, num_α, num_β)
    Σ = Array{Any}(undef, num_α, num_β, num_Ω)

    for i in eachindex(α)
        Mₖ[i] = pustrip(frohlich_coupling(α[i], ω * ω0_pu, mb; dims = dims))
        for j in eachindex(β)
            S(v, w) = frohlich_S(v, w, Mₖ[i], τ -> β[j] == Inf ? phonon_propagator(τ, ω) : phonon_propagator(τ, ω, β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω : polaron_propagator(τ, v, w, β[j]) * ω; dims = dims, limits = [0, sqrt(β[j] / 2)])
            v_guess = v_guesses == false ? α[i] < 7 ? 3 + α[i] / 4 + 1 / β[j] : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
            w_guess = w_guesses == false ? 2 + tanh((6 - α[i]) / 3) + 1 / β[j] : w_guesses
            v[i,j], w[i,j], E[i,j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[i,j,k] = isinf(β[j]) && Ω[k] == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ[i], t -> β[j] == Inf ? phonon_propagator(t, ω) : phonon_propagator(t, ω, β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,j], w[i,j]) * ω : polaron_propagator(t, v[i,j], w[i,j], β[j]) * ω; dims = dims) * ω
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Frohlich(ω * ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω .* ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::AbstractArray, Ω::Number; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)

    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_ω, num_β = length(α), length(ω), length(β)
    if verbose N, n = num_α * num_ω * num_β, 1 end

    Mₖ = Matrix{Any}(undef, num_α, num_ω)
    v = Array{Any}(undef, num_α, num_ω, num_β)
    w = Array{Any}(undef, num_α, num_ω, num_β)
    E = Array{Any}(undef, num_α, num_ω, num_β)
    Σ = Array{Any}(undef, num_α, num_ω, num_β)

    for i in eachindex(α), p in eachindex(ω)
        Mₖ[i,p] = pustrip(frohlich_coupling(α[i], ω[p] * ω0_pu, mb; dims = dims))
        for j in eachindex(β)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β]") end
            S(v, w) = frohlich_S(v, w, Mₖ[i,p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β[j]) * ω[p]; dims = dims, limits = [0, sqrt(β[j] / 2)])
            v_guess = v_guesses == false ? α[i] < 7 ? 3 + α[i] / 4 + 1 / β[j] : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
            w_guess = w_guesses == false ? 2 + tanh((6 - α[i]) / 3) + 1 / β[j] : w_guesses
            v[i,p,j], w[i,p,j], E[i,p,j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
            Σ[i,p,j] = isinf(β[j]) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω, Mₖ[i,p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,p,j], w[i,p,j]) * ω[p] : polaron_propagator(t, v[i,p,j], w[i,p,j], β[j]) * ω[p]; dims = dims) * ω[p]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::Number, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)

    ω = pustrip.(ω)
    β = pustrip(β * ħ_pu / E0_pu) 
    Ω = pustrip.(Ω)

    num_α, num_ω, num_Ω = length(α), length(ω), length(Ω)
    if verbose N, n = num_α * num_ω * num_Ω, 1 end

    Mₖ = Matrix{Any}(undef, num_α, num_ω)
    v = Matrix{Any}(undef, num_α, num_ω)
    w = Matrix{Any}(undef, num_α, num_ω)
    E = Matrix{Any}(undef, num_α, num_ω)
    Σ = Array{Any}(undef, num_α, num_ω, num_Ω)

    for i in eachindex(α), p in eachindex(ω)
        Mₖ[i,p] = pustrip(frohlich_coupling(α[i], ω[p] * ω0_pu, mb; dims = dims))
        S(v, w) = frohlich_S(v, w, Mₖ[i,p], τ -> β == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β), (τ, v, w) -> β == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β) * ω[p]; dims = dims, limits = [0, sqrt(β / 2)])
        v_guess = v_guesses == false ? α[i] < 7 ? 3 + α[i] / 4 + 1 / β : 4 * α[i]^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β : v_guesses
        w_guess = w_guesses == false ? 2 + tanh((6 - α[i]) / 3) + 1 / β : w_guesses
        v[i,p], w[i,p], E[i,p] = variation((v, w) -> β == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β) * dims / 3 - (S(v, w) - S₀(v, w, β) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
        for k in eachindex(Ω)
            if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | Ω = $(Ω[k]) [$k/$num_Ω]") end
            Σ[i,p,k] = isinf(β) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ[i,p], t -> β == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β), t -> β == Inf ? polaron_propagator(t, v[i,p], w[i,p]) * ω[p] : polaron_propagator(t, v[i,p], w[i,p], β) * ω[p]; dims = dims) * ω[p]
            if verbose n += 1; print("\e[1F") end
        end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β / E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
end

function frohlich(α::AbstractArray, ω::AbstractArray, β::AbstractArray, Ω::AbstractArray; mb = 1m0_pu, dims = 3, v_guesses = false, w_guesses = false, upper = Inf, lower = eps(), verbose = false)

    ω = pustrip.(ω)
    β = pustrip.(β .* ħ_pu / E0_pu) 
    Ω = pustrip(Ω)

    num_α, num_ω, num_β, num_Ω = length(α), length(ω), length(β), length(Ω)
    if verbose N, n = num_α * num_ω * num_β * num_Ω, 1 end

    Mₖ = Matrix{Any}(undef, num_α, num_ω)
    v = Array{Any}(undef, num_α, num_ω, num_β)
    w = Array{Any}(undef, num_α, num_ω, num_β)
    E = Array{Any}(undef, num_α, num_ω, num_β)
    Σ = Array{Any}(undef, num_α, num_ω, num_β, num_Ω)

    for i in eachindex(α), p in eachindex(ω)
        Mₖ[i,p] = pustrip(frohlich_coupling(α[i], ω[p] * ω0_pu, mb; dims = dims))
        for j in eachindex(β)
            S(v, w) = frohlich_S(v, w, Mₖ[i,p], τ -> β[j] == Inf ? phonon_propagator(τ, ω[p]) : phonon_propagator(τ, ω[p], β[j]), (τ, v, w) -> β[j] == Inf ? polaron_propagator(τ, v, w) * ω[p] : polaron_propagator(τ, v, w, β[j]) * ω[p]; dims = dims, limits = [0, sqrt(β[j] / 2)])
            v_guess = v_guesses == false ? α < 7 ? 3 + α / 4 + 1 / β[j] : 4 * α^2 / 9π - 3/2 * (2 * log(2) + 0.5772) - 3/4 + 1 / β[j] : v_guesses
            w_guess = w_guesses == false ? 2 + tanh((6 - α) / 3) + 1 / β[j] : w_guesses
            v[i,p,j], w[i,p,j], E[i,p,j] = variation((v, w) -> β[j] == Inf ? E₀(v, w) * dims / 3 - (S(v, w) - S₀(v, w) * dims / 3) : E₀(v, w, β[j]) * dims / 3 - (S(v, w) - S₀(v, w, β[j]) * dims / 3), v_guess, w_guess; upper = upper, lower = lower)
            for k in eachindex(Ω)
                if verbose println("\e[K[$n/$N ($(round(n/N*100, digits=1)) %)] | α = $(α[i]) [$i/$num_α] | ω = $(ω[p]) [$p/$num_ω] | β = $(β[j]) [$j/$num_β] | Ω = $(Ω[k]) [$k/$num_Ω]") end
                Σ[i,p,j,k] = isinf(β[j]) && Ω == eps() ? zero(Complex) : frohlich_memory(Ω[k], Mₖ[i,p], t -> β[j] == Inf ? phonon_propagator(t, ω[p]) : phonon_propagator(t, ω[p], β[j]), t -> β[j] == Inf ? polaron_propagator(t, v[i,p,j], w[i,p,j]) * ω[p] : polaron_propagator(t, v[i,p,j], w[i,p,j], β[j]) * ω[p]; dims = dims) * ω[p]
                if verbose n += 1; print("\e[1F") end
            end
        end
    end
    return Frohlich(ω .* ω0_pu, Mₖ .* E0_pu, α, E .* E0_pu, v .* ω0_pu, w .* ω0_pu, β ./ E0_pu, Ω * ω0_pu, Σ .* ω0_pu)
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

function optimal_frohlich(αωβΩ...; rtol = 1e-5, upper = Inf, lower = eps(), verbose = false, kwargs...) 
    if verbose println("Estimating n = 1 Fictitious Particle...") end
    args = [x[1] for x in αωβΩ]
    f = frohlich(args...; v_guesses = [3], w_guesses = [2], lower = lower, upper = upper, kwargs...)
    if verbose println("Energy: $(pustrip.(f.E)) | v parameters: $(pustrip.(f.v)) | w parameters: $(pustrip.(f.w))") end
    energy = f.E
    err = 1
    v_guesses, w_guesses = pustrip.(f.v), pustrip.(f.w)
    new_lower, new_upper = lower, upper
    if verbose n = 2 end
    while err > rtol
        if verbose println("Estimating n = $n Fictitious Particles...") end
        new_upper = pustrip.(vcat(repeat([f.w...], inner=2), [maximum(f.v) * 3, maximum(f.v) * 3]))
        v_guesses = pustrip.(vcat(f.v .* 0.5, maximum(f.v) * 2))
        w_guesses = pustrip.(vcat(f.w .* 0.5, maximum(f.w) * 2))
        new_lower = repeat([lower], length(v_guesses) + length(w_guesses))
        println(new_upper)
        println(new_lower)
        f = frohlich(args...; v_guesses = v_guesses, w_guesses = w_guesses, lower = new_lower, upper = new_upper, kwargs...)
        if verbose println("Energy: $(pustrip.(f.E)) | v parameters: $(pustrip.(f.v)) | w parameters: $(pustrip.(f.w))") end
        new_energy = f.E
        err = (new_energy - energy) / energy
        if verbose println("Error: $(maximum(err))") end
        energy = new_energy
        if verbose n += 1 end
    end
    if verbose 
        println("Optimal number of fictitious particles: $(n - 1)")
        println("Relative error: $(err)")
        println("Evaluating polaron properties...")
    end
    return f = frohlich(αωβΩ...; v_guesses = v_guesses, w_guesses = w_guesses, lower = new_lower, upper = new_upper, verbose = verbose, kwargs...)
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
    integral, _ = quadgk(τ -> 2 * τ * phonon_propagator(τ^2) / sqrt(polaron_propagator(τ^2, v, w)), limits...)
    return norm(coupling)^2 * ball_surface(dims) / (2π)^dims * sqrt(π / 2) * integral
end

function frohlich_memory(Ω, coupling, phonon_propagator, polaron_propagator; dims = 3)
    integral, _ = quadgk(t -> 2 * t * (1 - exp(im * Ω * t^2)) / Ω * imag(phonon_propagator(im * t^2) / polaron_propagator(im * t^2)^(3/2)), 0, Inf)
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
