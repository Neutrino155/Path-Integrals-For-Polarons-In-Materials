# Material.jl

abstract type Material end

struct FrohlichMaterial <: Material
    formula::String
    ϵ_total
    ϵ_optic
    ϵ_ionic
    ω_LO
    mb
    V
    function FrohlichMaterial(x...)
        new(x...)
    end
end

frohlichmaterial() = FrohlichMaterial([], [], [], [], [], [], [])

# frohlichmaterial(x...) = FrohlichMaterial(x...)

frohlichmaterial(formula, ϵ_total::Number, ϵ_optic::Number, ω_LO::Number, mb::Number, V::Number) = FrohlichMaterial(formula, ϵ_total, ϵ_optic, ϵ_total - ϵ_optic, ω_LO, mb, V)

frohlichmaterial(formula, ϵ_total::Number, ϵ_optic::Number, ϵ_ionic::AbstractVector, ω_LO::AbstractVector, mb::Number, V::Number) = FrohlichMaterial(formula, ϵ_total, ϵ_optic, ϵ_ionic, ω_LO, mb, V)

function save_frohlich_material(data::FrohlichMaterial, prefix)

    println("Saving Frohlich material data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "formula", data.formula,
        "ϵ_total", pustrip.(data.ϵ_total),
        "ϵ_optic", pustrip.(data.ϵ_optic),
        "ϵ_ionic", pustrip.(data.ϵ_ionic),
        "ω_LO", pustrip.(data.ω_LO),
        "mb", pustrip.(data.mb),
        "V", pustrip.(data.V)
        )

    println("... Frohlich material data saved.")
end

function load_frohlich_material(frohlich_file_path)

    println("Loading Frohlich material data from $frohlich_file_path ...")

    data = JLD.load("$frohlich_file_path")

    frohlich = FrohlichMaterial(
        data["formula"],
        data["ϵ_total"] .* punit(u"ϵ0"),
        data["ϵ_optic"] .* punit(u"ϵ0"),
        data["ϵ_ionic"] .* punit(u"ϵ0"),
        data["ω_LO"] .* ω0_pu,
        data["mb"] .* m0_pu,
        data["V"] .* r0_pu^3
    )
    println("... Frohlich material data loaded.")

    return frohlich
end

struct HolsteinMaterial <: Material
    formula             # Chemical formula of material
    J                   # Transfer integral
    ω_LO                # Phonon frequency
    g                   # El-ph coupling energy
    lattice_constant    # Lattice constant
    d                   # dimensionality
    function HolsteinMaterial(x...)
        new(x...)
    end
end

# holsteinmaterial = HolsteinMaterial([], [], [], [], [])

function save_holstein_material(data::HolsteinMaterial, prefix)

    println("Saving Holstein material data to $prefix.jld ...")

    JLD.save("$prefix.jld",
        "formula", data.formula,
        "J", pustrip.(data.J),
        "ω", pustrip.(data.ω_LO),
        "g", pustrip.(data.g),
        "a", pustrip.(data.lattice_constant),
        "d", pustrip.(data.d)
        )

    println("... Holstein material data saved.")
end

function load_holstein_material(holstein_file_path)

    println("Loading Holstein material data from $holstein_file_path ...")

    data = JLD.load("$holstein_file_path")

    holstein = HolsteinMaterial(
        data["formula"],
        data["J"] .* E0_pu,
        data["ω"] .* ω0_pu,
        data["g"] .* E0_pu,
        data["a"] .* r0_pu,
        data["d"] 
    )
    println("... Holstein material data loaded.")

    return holstein
end