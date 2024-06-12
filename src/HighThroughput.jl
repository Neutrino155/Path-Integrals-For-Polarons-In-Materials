# HighThroughput.jl

# Standard and feynman data

standard_data_conduction = readdlm("data/LiegeDataset/Results/StandardFeynmanFrohlich/valence/standard_feynman_data_valence")

for material in eachrow(standard_data_conduction)[2:end]
    println(material[1])
    volume = parse(Float64, split(split(JSON.parsefile("data/LiegeDataset/Repository/phonon/$(material[1]).json")["metadata"]["structure"], "\n")[13], "   ")[2])
    f = frohlichmaterial(material[2],  material[3] * u"ϵ0", material[4] * u"ϵ0", material[5] * u"meV", material[6]^2 * u"me", volume * u"Å^3")
    save_frohlich_material(f, "data/frohlich/materials/feynman/valence/$(material[1])-$(material[2])")
end

# General Materials

generalmaterials = readdir("data/LiegeDataset/Results/GeneralizedFrohlich/conduction/")[1:1395]

for file in generalmaterials

    mp_id = join(split(file, "-")[2:3], "-")
    material = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/conduction/$file")[2,5:7]
    formula = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/conduction/$file")[2,2] 
    
    if !iszero(material[2]) && !iszero(material[3])

        volume = parse(Float64, split(split(JSON.parsefile("data/LiegeDataset/Repository/phonon/$mp_id.json")["metadata"]["structure"], "\n")[13], "   ")[2]) * u"Å^3"
        ϵ_optic = det(hcat(JSON.parsefile("data/LiegeDataset/Repository/phonon/$mp_id.json")["dielectric"]["eps_electronic"]...))^(1/3) * u"ϵ0"
        mb = material[1]^2 * u"me"

        α = material[2]
        E = material[3] * u"meV"
        ω = abs(E / α) / ħ_pu
        rp = sqrt(ħ_pu / (2 * mb * ω))
        ϵ_total = (1/ϵ_optic - 8π * ħ_pu * ω * α * rp / e_pu^2)^(-1)

        f_single = frohlichmaterial(formula, ϵ_total, ϵ_optic, ω, mb, volume)
        save_frohlich_material(f_single, "data/frohlich/materials/general/conduction/$mp_id-$formula-single")
    end
end

generalmaterials = readdir("data/LiegeDataset/Results/GeneralizedFrohlich/conduction/")[1396:end]

for file in generalmaterials

    mp_id = join(split(file, "-")[1:2], "-")
    material = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/conduction/$file")[2:end, 4:5]

    if isfile("data/LiegeDataset/Results/GeneralizedFrohlich/conduction/gFr-$(mp_id)-$(split(file, "mode-")[2])")
        material_single = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/conduction/gFr-$(mp_id)-$(split(file, "mode-")[2])")[2, :]
        formula = material_single[2]

        multimode = material[[all(row .> 0.00) for row in eachrow(abs.(material))], :]

        volume = parse(Float64, split(split(JSON.parsefile("data/LiegeDataset/Repository/phonon/$mp_id.json")["metadata"]["structure"], "\n")[13], "   ")[2]) * u"Å^3"
        ϵ_optic = det(hcat(JSON.parsefile("data/LiegeDataset/Repository/phonon/$mp_id.json")["dielectric"]["eps_electronic"]...))^(1/3) * u"ϵ0"
        mb = material_single[5]^2 * u"me"

        if size(multimode)[1] > 1

            α = multimode[:, 1]
            E = multimode[:, 2] .* u"meV"
            ω = abs.(E ./ α) ./ ħ_pu
            rp = sqrt.(ħ_pu ./ (2 * mb .* ω))
            ϵ_total = abs(1/ϵ_optic - sum(8π .* ħ_pu .* ω .* α .* rp ./ e_pu^2))^(-1)
            ϵ_ionic = 8π .* α .* ħ_pu .* ω .* rp .* ϵ_optic .* ϵ_total ./ e_pu^2

            f_multi = frohlichmaterial(formula, ϵ_total, ϵ_optic, ϵ_ionic, ω, mb, volume)
            save_frohlich_material(f_multi, "data/frohlich/materials/general/conduction/$mp_id-$formula-multi")
        end
    end
end

# Generate Data

frohlichmaterials = readdir("data/frohlich/materials/general/conduction/")

for material in frohlichmaterials

    file_name = split(material, ".")[1]

    m = load_frohlich_material("data/frohlich/materials/general/conduction/$material")

    # f_gs = frohlich(m)
    # save_frohlich(f_gs, "data/frohlich/variational/general/valence/groundstate/$file_name-gs")

    f_rt = frohlich(m, 300u"K")
    save_frohlich(f_rt, "data/frohlich/variational/general/conduction/roomtemp/$file_name-rt")
end


# Collate data

materials = readdir("data/frohlich/materials/feynman/valence")
liege = readdlm("data/LiegeDataset/Results/StandardFeynmanFrohlich/valence/standard_feynman_data_valence")

headers = Vector{Any}(["mp_id", "formula", "eps_total(ϵ0)", "eps_optic(ϵ0)", "eps_ionic(ϵ0)", "omega_LO(THz2π)", "band_mass(me)", "volume(Å^3)", "alpha", "v_gs", "w_gs", "v_rt", "w_rt", "energy_gs(meV)", "energy_rt(meV)", "ZPR(meV)", "imag_memory", "mobility(cm^2/V/s)"])

data_all = [headers]

for material in materials
    file_name = split(material, ".")[1]
    mp_id = join(split(file_name, "-")[1:2], "-")
    formula = split(file_name, "-")[3]
    if isfile("data/frohlich/variational/feynman/valence/roomtemp/$(file_name)-rt.jld")
        m = load_frohlich_material("data/frohlich/materials/feynman/valence/$(material)")
        f_gs = load_frohlich("data/frohlich/variational/feynman/valence/groundstate/$(file_name)-gs.jld")
        f_rt = load_frohlich("data/frohlich/variational/feynman/valence/roomtemp/$(file_name)-rt.jld")
        data = [mp_id, formula, m.ϵ_total / u"ϵ0" |> NoUnits, m.ϵ_optic / u"ϵ0" |> NoUnits, m.ϵ_ionic / u"ϵ0" |> NoUnits, m.ω_LO / ω0_pu, m.mb / m0_pu, m.V / u"Å^3" |> NoUnits, f_gs.α, f_gs.v / f_gs.ω, f_gs.w / f_gs.ω, f_rt.v / f_rt.ω, f_rt.w / f_rt.ω, ustrip(f_gs.E |> u"meV"), ustrip(f_rt.E |> u"meV"), liege[liege[:,1] .== mp_id,end][1], imag.(f_rt.Σ) / ω0_pu, ustrip(1 ./ imag.(f_rt.Σ) * e_pu / m.mb |> u"cm^2/V/s")]
        push!(data_all, data)
    end
end

writedlm("frohlich-materials-feynman-valence.txt", data_all)


# general collate data

materials = readdir("data/frohlich/materials/general/valence")

headers = Vector{Any}(["mp_id", "formula", "mode", "eps_total(ϵ0)", "eps_optic(ϵ0)", "eps_ionic(ϵ0)", "band_mass(me)", "volume(Å^3)", "alpha", "v_gs", "w_gs", "v_rt", "w_rt", "energy_gs(meV)", "energy_rt(meV)", "AHC_alpha", "AHC_ZPR(meV)", "imag_memory", "mobility(cm^2/V/s)"])

data_all = [headers]

for material in materials
    file_name = split(material, ".")[1]
    mp_id = join(split(file_name, "-")[1:2], "-")
    mode = split(file_name, "-")[4]
    formula = split(file_name, "-")[3]

    m = load_frohlich_material("data/frohlich/materials/general/valence/$file_name.jld")
    f_gs = load_frohlich("data/frohlich/variational/general/valence/groundstate/$file_name-gs.jld")
    f_rt = load_frohlich("data/frohlich/variational/general/valence/roomtemp/$file_name-rt.jld")

    if mode == "single"
        liege = isfile("data/LiegeDataset/Results/GeneralizedFrohlich/valence/gFr-$(mp_id)-valence.dat") ? readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/valence/gFr-$(mp_id)-valence.dat")[2:end, :] : readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/valence/gFr-$(mp_id)-valence-1.dat")[2:end, :]
        liege_alpha = liege[6]
        liege_alpha_sum = sum(liege_alpha)
        liege_ZPR = liege[7]
    else
        liege = isfile("data/LiegeDataset/Results/GeneralizedFrohlich/valence/$(mp_id)-data-per-mode-valence.dat") ? readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/valence/$(mp_id)-data-per-mode-valence.dat")[2:end, :] : readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/valence/$(mp_id)-data-per-mode-valence-1.dat")[2:end, :]
        liege_alpha = liege[:, 4]
        liege_alpha_sum = sum(liege_alpha)
        liege_ZPR = sum(liege[:, 5])
    end

    data = [mp_id, formula, mode, m.ϵ_total ./ u"ϵ0" .|> NoUnits, m.ϵ_optic ./ u"ϵ0" .|> NoUnits, sum(m.ϵ_ionic) ./ u"ϵ0" .|> NoUnits, m.mb ./ m0_pu, m.V ./ u"Å^3" .|> NoUnits, sum(f_gs.α), f_gs.v ./ ω0_pu, f_gs.w ./ ω0_pu, f_rt.v ./ ω0_pu, f_rt.w ./ ω0_pu, ustrip.(f_gs.E .|> u"meV"), ustrip.(f_rt.E .|> u"meV"), liege_alpha_sum, liege_ZPR, imag.(f_rt.Σ) ./ ω0_pu, ustrip(1 ./ imag.(f_rt.Σ) .* e_pu ./ m.mb .|> u"cm^2/V/s")]

    push!(data_all, data)
end

writedlm("frohlich-materials-general-valence.txt", data_all)

