### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 9d01a510-3ec7-11ef-0f8b-7be0ed93f14d
using Pkg, JSON, DelimitedFiles, Unitful, LinearAlgebra

# ╔═╡ b557b342-c459-4dfc-9ffe-aa58cf46286c
Pkg.activate(".")

# ╔═╡ 6a84ff3e-40b9-46d7-9e94-33a4db00d5ac
using PolaronThesis

# ╔═╡ f3d30776-d86d-4f59-b200-a60cbcfb0f2b
element_groups = Dict(
    "H" => 1, "He" => 18, "Li" => 1, "Be" => 2, "B" => 13, "C" => 14, "N" => 15, "O" => 16, "F" => 17, "Ne" => 18,
    "Na" => 1, "Mg" => 2, "Al" => 13, "Si" => 14, "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 1, "Ca" => 2,
    "Sc" => 3, "Ti" => 4, "V" => 5, "Cr" => 6, "Mn" => 7, "Fe" => 8, "Co" => 9, "Ni" => 10, "Cu" => 11, "Zn" => 12,
    "Ga" => 13, "Ge" => 14, "As" => 15, "Se" => 16, "Br" => 17, "Kr" => 18, "Rb" => 1, "Sr" => 2, "Y" => 3, "Zr" => 4,
    "Nb" => 5, "Mo" => 6, "Tc" => 7, "Ru" => 8, "Rh" => 9, "Pd" => 10, "Ag" => 11, "Cd" => 12, "In" => 13, "Sn" => 14,
    "Sb" => 15, "Te" => 16, "I" => 17, "Xe" => 18, "Cs" => 1, "Ba" => 2, "La" => 3, "Ce" => 4, "Pr" => 5, "Nd" => 6,
    "Pm" => 7, "Sm" => 8, "Eu" => 9, "Gd" => 10, "Tb" => 11, "Dy" => 12, "Ho" => 13, "Er" => 14, "Tm" => 15, "Yb" => 16,
    "Lu" => 17, "Hf" => 4, "Ta" => 5, "W" => 6, "Re" => 7, "Os" => 8, "Ir" => 9, "Pt" => 10, "Au" => 11, "Hg" => 12,
    "Tl" => 13, "Pb" => 14, "Bi" => 15, "Po" => 16, "At" => 17, "Rn" => 18, "Fr" => 1, "Ra" => 2, "Ac" => 3, "Th" => 4,
    "Pa" => 5, "U" => 6, "Np" => 7, "Pu" => 8, "Am" => 9, "Cm" => 10, "Bk" => 11, "Cf" => 12, "Es" => 13, "Fm" => 14,
    "Md" => 15, "No" => 16, "Lr" => 17
)

# ╔═╡ d18b8069-2dd5-402e-bd20-bebfee3bc0b7
function isnumeric(str)
    # Check if a string represents a number
    try
        parse(Int, str)
        return true
    catch
        return false
    end
end

# ╔═╡ 8ca72347-d4e4-4924-8f3c-740825d55b67
function split_formula(formula)
    # Split the formula at every capital letter
    parts = split(formula, r"(?=[A-Z])")
    
    # Separate elements and numbers
    result = []
    for part in parts
        # Split numbers and elements within each part
        sub_parts = split(part, r"(?=\d)")
        for sub_part in sub_parts
            if isnumeric(sub_part)
                push!(result, parse(Int, sub_part))
            else
                push!(result, sub_part)
            end
        end
    end
    return result
end

# ╔═╡ 8ea70512-c872-427a-8367-865dee88da06
# Function to get groups from a chemical formula
function get_element_groups(formula)
    elements = split_formula(formula)
    groups = []
    for elem in elements
        if elem isa SubString
            if haskey(element_groups, elem)
                push!(groups, element_groups[elem])
            else
                println("Warning: Element $elem not found in the periodic table.")
            end
        end
    end
    return groups
end

# ╔═╡ cb752228-e2fe-478f-a495-ed16ba97be45
# Standard and Feynman Materials

# ╔═╡ 2b339d71-91f5-4058-8973-c99d326c918e
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"]
	standard_feynman_data = readdlm("data/LiegeDataset/Results/StandardFeynmanFrohlich/$band/standard_feynman_data_$band")	
	for material in eachrow(standard_feynman_data)[2:end]
	    volume = parse(Float64, split(split(JSON.parsefile("data/LiegeDataset/Repository/phonon/$(material[1]).json")["metadata"]["structure"], "\n")[13], "   ")[2])
	    f = frohlichmaterial(material[2],  material[3] * u"ϵ0", material[4] * u"ϵ0", material[5] * u"meV", material[6]^2 * u"me", volume * u"Å^3")
	    save_frohlich_material(f, "data/frohlich/materials/feynman/$band/$(material[1])-$(material[2])")
	end
end
  ╠═╡ =#

# ╔═╡ 0e6c408d-5a4f-4577-8482-91681761d5c6
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"]
	standard_data = readdlm("data/LiegeDataset/Results/StandardFrohlich/$band/standard_froelich_data_$band")	
	for material in eachrow(standard_data)[2:end]
	    volume = parse(Float64, split(split(JSON.parsefile("data/LiegeDataset/Repository/phonon/$(material[1]).json")["metadata"]["structure"], "\n")[13], "   ")[2])
	    f = frohlichmaterial(material[2],  material[3] * u"ϵ0", material[4] * u"ϵ0", material[5] * u"meV", material[6]^2 * u"me", volume * u"Å^3")
	    save_frohlich_material(f, "data/frohlich/materials/standard/$band/$(material[1])-$(material[2])")
	end
end
  ╠═╡ =#

# ╔═╡ 6568d544-c9bf-4722-a5cf-3aca5e56e9d2
# General Materials

# ╔═╡ adb00fa6-9f7f-4f69-9134-4507de6dfe36
# Single mode materials

# ╔═╡ 52d2f68e-1e36-4f8c-8df3-657d364953e4
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"]
	
	singlemode_materials = readdir("data/LiegeDataset/Results/GeneralizedFrohlich/$band/")[1:1394]

	for file in singlemode_materials

	    mp_id = join(split(file, "-")[2:3], "-")
	    material = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/$file")[2,5:7]
	    formula = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/$file")[2,2] 
	    
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
	        save_frohlich_material(f_single, "data/frohlich/materials/general/$band/$mp_id-$formula-single")
	    end
	end
end
  ╠═╡ =#

# ╔═╡ c19253d3-da29-4a06-852f-742bd1d1f13e
# Multiple mode materials

# ╔═╡ 81ff1db1-b8b6-4462-a00f-88550669f8ce
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"]

	multimode_materials = readdir("data/LiegeDataset/Results/GeneralizedFrohlich/$band/")[1395:end]

	for file in multimode_materials

	    mp_id = join(split(file, "-")[1:2], "-")
	    material = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/$file")[2:end, 4:5]
	
	    if isfile("data/LiegeDataset/Results/GeneralizedFrohlich/$band/gFr-$(mp_id)-$(split(file, "mode-")[2])")
			
	        material_single = readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/gFr-$(mp_id)-$(split(file, "mode-")[2])")[2, :]
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
	            save_frohlich_material(f_multi, "data/frohlich/materials/general/$band/$mp_id-$formula-multi")
	        end
	    end
	end
end
  ╠═╡ =#

# ╔═╡ b8c34cb1-06f5-4bd1-a4c3-4d90ea235fb9
# Generate Frohlich polaron data

# ╔═╡ 8dfe9b3f-4c7f-4d46-9301-7a06e7297bf6
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"], theory in ["standard", "feynman", "general"]
	frohlich_materials = readdir("data/frohlich/materials/$theory/$band/")

	for material in frohlich_materials
	
	    file_name = split(material, ".")[1]
	
	    m = load_frohlich_material("data/frohlich/materials/$theory/$band/$material")
	
	    f_gs = frohlich(m)
	    save_frohlich(f_gs, "data/frohlich/variational/$theory/$band/groundstate/$file_name-gs")
	
	    f_rt = frohlich(m, 300u"K", eps())
	    save_frohlich(f_rt, "data/frohlich/variational/$theory/$band/roomtemp/$file_name-rt")
	end
end
  ╠═╡ =#

# ╔═╡ 6bc9fd00-a625-40cd-896e-a908c993305d
# Collate data into single files

# ╔═╡ 42ecd229-3cfc-4060-8476-ab6e1a29c989
# Standard and Feynman

# ╔═╡ 0c4cc6b0-185b-45c8-9014-4f939e3c7c54
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"], theory in ["standard", "feynman"]
	materials = readdir("data/frohlich/materials/$theory/$band")
	liege = theory == "standard" ? readdlm("data/LiegeDataset/Results/StandardFrohlich/$band/standard_froelich_data_$band") : readdlm("data/LiegeDataset/Results/StandardFeynmanFrohlich/$band/standard_feynman_data_$band")
	
	headers = Vector{Any}(["mp_id", "formula", "spgroup", "groups", "eps_total(ϵ0)", "eps_optic(ϵ0)", "eps_ionic(ϵ0)", "omega_LO(THz2π)", "band_mass(me)", "volume(Å^3)", "alpha", "v_gs", "w_gs", "v_rt", "w_rt", "energy_gs(meV)", "energy_rt(meV)", "ZPR(meV)", "imag_memory", "mobility(cm^2/V/s)"])
	
	data_all = [headers]
	
	for material in materials
	    file_name = split(material, ".")[1]
	    mp_id = join(split(file_name, "-")[1:2], "-")
	    formula = split(file_name, "-")[3]
		groups = join(get_element_groups(formula), "-")
		abinit_data = readdlm("data/LiegeDataset/Repository/abinit_effmasses/$(mp_id)_output")
	    spgroup = abinit_data[findfirst(x -> x == "spgroup", abinit_data[:,1]),2]
		
	    if isfile("data/frohlich/variational/$theory/$band/roomtemp/$(file_name)-rt.jld")
	
			m = load_frohlich_material("data/frohlich/materials/$theory/$band/$(material)")
			f_gs = load_frohlich("data/frohlich/variational/$theory/$band/groundstate/$(file_name)-gs.jld")
			f_rt = load_frohlich("data/frohlich/variational/$theory/$band/roomtemp/$(file_name)-rt.jld")
	
			data = [mp_id, formula, spgroup, groups, m.ϵ_total / u"ϵ0" |> NoUnits, m.ϵ_optic / u"ϵ0" |> NoUnits, m.ϵ_ionic / u"ϵ0" |> NoUnits, m.ω_LO / ω0_pu, m.mb / m0_pu, m.V / u"Å^3" |> NoUnits, f_gs.α, f_gs.v / f_gs.ω, f_gs.w / f_gs.ω, f_rt.v / f_rt.ω, f_rt.w / f_rt.ω, ustrip(f_gs.E |> u"meV"), ustrip(f_rt.E |> u"meV"), liege[liege[:,1] .== mp_id,end][1], imag.(f_rt.Σ) / ω0_pu, ustrip(1 ./ imag.(f_rt.Σ) * e_pu / m.mb |> u"cm^2/V/s")]
				
			push!(data_all, data)
			
	    end
	end
	
	writedlm("data/frohlich/variational/$theory/frohlich-materials-$theory-$band.txt", data_all)
end
  ╠═╡ =#

# ╔═╡ e0d2929d-5b42-467b-aef8-5bfa0674be38
# General

# ╔═╡ c729301e-000b-49d0-91df-22f57331cdc4
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"]
	
	materials = readdir("data/frohlich/materials/general/$band")

	headers = Vector{Any}(["mp_id", "formula", "spgroup", "groups", "mode", "eps_total(ϵ0)", "eps_optic(ϵ0)", "eps_ionic(ϵ0)", "band_mass(me)", "volume(Å^3)", "alpha", "v_gs", "w_gs", "v_rt", "w_rt", "energy_gs(meV)", "energy_rt(meV)", "AHC_alpha", "AHC_ZPR(meV)", "imag_memory", "mobility(cm^2/V/s)"])
	
	data_all = [headers]
	
	for material in materials
	    file_name = split(material, ".")[1]
	    mp_id = join(split(file_name, "-")[1:2], "-")
	    mode = split(file_name, "-")[4]
	    formula = split(file_name, "-")[3]
		groups = join(get_element_groups(formula), "-")
		abinit_data = readdlm("data/LiegeDataset/Repository/abinit_effmasses/$(mp_id)_output")
	    spgroup = abinit_data[findfirst(x -> x == "spgroup", abinit_data[:,1]),2]
	
	    m = load_frohlich_material("data/frohlich/materials/general/$band/$file_name.jld")
	    f_gs = load_frohlich("data/frohlich/variational/general/$band/groundstate/$file_name-gs.jld")
	    f_rt = load_frohlich("data/frohlich/variational/general/$band/roomtemp/$file_name-rt.jld")
	
	    if mode == "single"
	        liege = isfile("data/LiegeDataset/Results/GeneralizedFrohlich/$band/gFr-$(mp_id)-$band.dat") ? readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/gFr-$(mp_id)-$band.dat")[2:end, :] : readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/gFr-$(mp_id)-$band-1.dat")[2:end, :]
	        liege_alpha = liege[6]
	        liege_alpha_sum = sum(liege_alpha)
	        liege_ZPR = liege[7]
	    else
	        liege = isfile("data/LiegeDataset/Results/GeneralizedFrohlich/$band/$(mp_id)-data-per-mode-$band.dat") ? readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/$(mp_id)-data-per-mode-$band.dat")[2:end, :] : readdlm("data/LiegeDataset/Results/GeneralizedFrohlich/$band/$(mp_id)-data-per-mode-$band-1.dat")[2:end, :]
	        liege_alpha = liege[:, 4]
	        liege_alpha_sum = sum(liege_alpha)
	        liege_ZPR = sum(liege[:, 5])
	    end
	
	    data = [mp_id, formula, spgroup, groups, mode, m.ϵ_total ./ u"ϵ0" .|> NoUnits, m.ϵ_optic ./ u"ϵ0" .|> NoUnits, sum(m.ϵ_ionic) ./ u"ϵ0" .|> NoUnits, m.mb ./ m0_pu, m.V ./ u"Å^3" .|> NoUnits, sum(f_gs.α), f_gs.v ./ ω0_pu, f_gs.w ./ ω0_pu, f_rt.v ./ ω0_pu, f_rt.w ./ ω0_pu, ustrip.(f_gs.E .|> u"meV"), ustrip.(f_rt.E .|> u"meV"), liege_alpha_sum, liege_ZPR, imag.(f_rt.Σ) ./ ω0_pu, ustrip(1 ./ imag.(f_rt.Σ) .* e_pu ./ m.mb .|> u"cm^2/V/s")]
	
	    push!(data_all, data)
	end
	
	writedlm("data/frohlich/variational/general/frohlich-materials-general-$band.txt", data_all)
end
  ╠═╡ =#

# ╔═╡ 880e931f-0970-4d54-adb9-88383709493c
# ╠═╡ disabled = true
#=╠═╡
for band in ["conduction", "valence"]
	
	materials = eachrow(readdlm("data/frohlich/variational/general/frohlich-materials-general-$band.txt"))[2:end]

	headers = Vector{Any}(["mp_id", "formula", "spgroup", "groups", "volume(Å^3)", "band_mass(me)", "eps_optic(ϵ0)", "eps_total_single(ϵ0)", "eps_total_multi(ϵ0)", "eps_ionic_single(ϵ0)", "eps_ionic_multi(ϵ0)", "AHC_alpha", "alpha", "v_gs_single", "v_gs_multi", "w_gs_single", "w_gs_multi", "v_rt_single", "v_rt_multi", "w_rt_single", "w_rt_multi", "AHC_ZPR(meV)", "energy_gs_single(meV)", "energy_gs_multi(meV)", "energy_rt_single(meV)", "energy_rt_multi(meV)", "imag_memory_single", "imag_memory_multi", "mobility_single(cm^2/V/s)", "mobility_multi(cm^2/V/s)"])
	
	data_all = [headers]
	
	for material in eachindex(materials)
		if materials[material][5] == "multi"
			mp_id = materials[material][1]
			formula = materials[material][2]
			spgroup = materials[material][3]
			groups = materials[material][4]
			volume = materials[material][10]
			bandmass = materials[material][9]
			optic = materials[material][7]
			AHCalpha = materials[material][18]
			alpha = materials[material][11]
			ZPR = materials[material][19]

			total_single = materials[material+1][6]
			total_multi = materials[material][6]
			ionic_single = materials[material+1][8]
			ionic_multi = materials[material][8]
			vgs_single = materials[material+1][12]
			vgs_multi = materials[material][12]
			wgs_single = materials[material+1][13]
			wgs_multi = materials[material][13]
			vrt_single = materials[material+1][14]
			vrt_multi = materials[material][14]
			wrt_single = materials[material+1][15]
			wrt_multi = materials[material][15]
			Egs_single = materials[material+1][16]
			Egs_multi = materials[material][16]
			Ert_single = materials[material+1][17]
			Ert_multi = materials[material][17]
			imag_single = materials[material+1][20]
			imag_multi = materials[material][20]
			mobility_single = materials[material+1][21]
			mobility_multi = materials[material][21]

			data = [mp_id, formula, spgroup, groups, volume, bandmass, optic, total_single, total_multi, ionic_single, ionic_multi, AHCalpha, alpha, vgs_single, vgs_multi, wgs_single, wgs_multi, vrt_single, vrt_multi, wrt_single, wrt_multi, ZPR, Egs_single, Egs_multi, Ert_single, Ert_multi, imag_single, imag_multi, mobility_single, mobility_multi]

			push!(data_all, data)
		end
	end
	
	writedlm("data/frohlich/variational/general/frohlich-materials-general-combined-$band.txt", data_all)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═9d01a510-3ec7-11ef-0f8b-7be0ed93f14d
# ╠═b557b342-c459-4dfc-9ffe-aa58cf46286c
# ╠═6a84ff3e-40b9-46d7-9e94-33a4db00d5ac
# ╟─f3d30776-d86d-4f59-b200-a60cbcfb0f2b
# ╟─8ca72347-d4e4-4924-8f3c-740825d55b67
# ╟─d18b8069-2dd5-402e-bd20-bebfee3bc0b7
# ╟─8ea70512-c872-427a-8367-865dee88da06
# ╠═cb752228-e2fe-478f-a495-ed16ba97be45
# ╟─2b339d71-91f5-4058-8973-c99d326c918e
# ╟─0e6c408d-5a4f-4577-8482-91681761d5c6
# ╠═6568d544-c9bf-4722-a5cf-3aca5e56e9d2
# ╠═adb00fa6-9f7f-4f69-9134-4507de6dfe36
# ╟─52d2f68e-1e36-4f8c-8df3-657d364953e4
# ╠═c19253d3-da29-4a06-852f-742bd1d1f13e
# ╟─81ff1db1-b8b6-4462-a00f-88550669f8ce
# ╠═b8c34cb1-06f5-4bd1-a4c3-4d90ea235fb9
# ╟─8dfe9b3f-4c7f-4d46-9301-7a06e7297bf6
# ╠═6bc9fd00-a625-40cd-896e-a908c993305d
# ╠═42ecd229-3cfc-4060-8476-ab6e1a29c989
# ╟─0c4cc6b0-185b-45c8-9014-4f939e3c7c54
# ╠═e0d2929d-5b42-467b-aef8-5bfa0674be38
# ╠═c729301e-000b-49d0-91df-22f57331cdc4
# ╠═880e931f-0970-4d54-adb9-88383709493c
