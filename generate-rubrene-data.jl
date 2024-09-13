### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 69067b24-4916-11ef-037f-0902e798e948
using Pkg, Unitful, Plots, DelimitedFiles

# ╔═╡ ad02f2ba-b5d6-44ba-afcd-279be120007a
Pkg.activate(".")

# ╔═╡ 50109bde-f4ee-4918-9c74-753cad69a2e6
using PolaronThesis

# ╔═╡ bca41c3b-7daa-4445-997e-b3fef2e67d56
md"""
# Polaron Property Functions
"""

# ╔═╡ 1f1e81eb-a571-46c2-8618-e252cf69d205
rubrene_memory(h) = h.Σ .|> u"ps^-1"

# ╔═╡ 17068899-f58f-4415-9057-06531b35ab66
md"""
# Rubrene Literature Data
"""

# ╔═╡ 762f3da9-af65-4b10-8210-6c4255ad9526
rubrene_phonon_frequencies = [
	57.8,
	59.6,
	89.0,
	107.3,
	139.1,
	639.1,
	1011.2,
	1344.7,
	1593.3
] .* u"c * cm^-1" .|> u"THz2π"

# ╔═╡ 82184ac9-eebd-4a23-8347-965206153a16
rubrene_holstein_coupling_energies = [
	-1.7,
	1.4,
	1.6,
	-0.14,
	-2.3,
	-7.5,
	-3.6,
	19.8,
	-42.0
] .* u"meV"

# ╔═╡ fe1b43ac-be1c-4a77-bee3-dc3d5b227e2e
rubrene_peierls_coupling_energies = [
	0.85,
	-0.83,
	-4.8,
	2.8,
	-3.7,
	1.0,
	-0.04,
	0.04,
	-0.12
] .* u"meV"

# ╔═╡ c46076de-52ee-48ed-ab68-e07eefc168f3
rubrene_holstein_binding_energy = sum(rubrene_holstein_coupling_energies.^2 ./ (ħ_pu .* rubrene_phonon_frequencies)) |> u"meV"

# ╔═╡ 662636a3-c9b4-4552-96ce-360d3e0fc93d
rubrene_peierls_binding_energy = sum(rubrene_peierls_coupling_energies.^2 ./ (ħ_pu .* rubrene_phonon_frequencies)) |> u"meV"

# ╔═╡ 2c9625e7-11e5-4d56-84e8-3d15ba0a6c6a
rubrene_holstein_phonon_frequency = sum(rubrene_holstein_coupling_energies .^2) / rubrene_holstein_binding_energy / ħ_pu |> u"cm^-1 * c"

# ╔═╡ 0ba581f6-1e4c-4960-8668-c00c1084bd64
rubrene_peierls_phonon_frequency = sum(rubrene_peierls_coupling_energies .^2) / rubrene_peierls_binding_energy / ħ_pu |> u"cm^-1 * c"

# ╔═╡ 5f8a7623-bcb4-4797-98a9-d26450d39409
rubrene_phonon_frequency = sum(rubrene_peierls_coupling_energies .^2 .+ rubrene_holstein_coupling_energies .^2) / (rubrene_peierls_binding_energy .+ rubrene_holstein_binding_energy) / ħ_pu |> u"cm^-1 * c"

# ╔═╡ 6298f629-457d-4881-82c4-2c356de25079
rubrene_holstein_coupling_energy = sqrt(rubrene_holstein_binding_energy * ħ_pu * rubrene_holstein_phonon_frequency) |> u"meV"

# ╔═╡ 78b458f9-8c59-4949-9555-64a1eda46709
rubrene_peierls_coupling_energy = sqrt(rubrene_peierls_binding_energy * ħ_pu * rubrene_peierls_phonon_frequency) |> u"meV"

# ╔═╡ 9060c1ad-5ce1-4bd5-a2af-45af9c04f0c3
rubrene_transfer_integral = 134u"meV"

# ╔═╡ d71db010-e118-4f22-99a3-8327e1b0e303
rubrene_lattice_constant = 0.7193u"nm"

# ╔═╡ a1d005ec-94cf-48f4-a667-ea9733497bab
rubrene_band_mass = puconvert(ħ_pu^2 / (2 * rubrene_transfer_integral * 
rubrene_lattice_constant^2))

# ╔═╡ e6f4c89c-32cb-419f-9a92-fd4837632b81
rubrene_mass(h) = h.v.^2 ./ h.w.^2 .* rubrene_band_mass

# ╔═╡ 6d58793a-126c-44d5-adc7-b4acfbdae12a
rubrene_spring(h) = (h.v.^2 .- h.w.^2) .* rubrene_band_mass .|> u"N/m"

# ╔═╡ 0135db00-7d9a-45ec-b9b9-bd4359a02089
rubrene_mobility(holstein) = abs.(1 ./ imag.(holstein.Σ)) * e_pu / rubrene_band_mass .|> u"cm^2/V/s"

# ╔═╡ 3851a3fe-6f3b-4099-b7a3-038dd1f5bf63
md"""
# Rubrene Material Data
"""

# ╔═╡ e74f64c8-bedc-461e-9cf0-407f8f1d4313
md"""
## Single Mode
"""

# ╔═╡ cbf984e6-602b-4753-9ecc-db2ddbe28a6b
#=╠═╡
function rubrene_radius(h)
	v, w, d, β = reduce_array(h.v), reduce_array(h.w), reduce_array(h.d), reduce_array(h.β)
	return d / 2 .* v ./ (v.^2 .- w.^2) .* coth.(pustrip.(β .* v ./ 2)) .* sqrt(ħ_pu * rubrene_single.ω_LO / rubrene_band_mass) .|> u"Å"
end
  ╠═╡ =#

# ╔═╡ 799e2865-eb50-49ca-8aaf-e1d648ae2b25
#=╠═╡
rubrene_conductivity(h) = -im ./ (h.Ω .+ h.Σ') .* e_pu^2 / rubrene_band_mass / (ħ_pu / rubrene_band_mass / rubrene_single.ω_LO) .|> u"μS"
  ╠═╡ =#

# ╔═╡ 39eda869-dfab-47e7-9c8d-2ca734eed33e
md"""
## Multiple Modes
"""

# ╔═╡ 5270ea61-c216-47e8-a9de-87635dabecec
md"""
# Rubrene | T = 1:400 K 
"""

# ╔═╡ c200c995-4eb8-4347-82fc-7dc170401df8
md"""
## Single Mode
"""

# ╔═╡ 8acd1416-755e-405f-852f-884a7976ca82
#=╠═╡
begin
data1 = Matrix{Any}(undef, length(rubrene_single_holstein_temp_1to400K.β) + 1, 10)
data1[1, :] = ["T(K)" "v(THz2π)" "w(THz2π)" "E(meV)" "M(mₑ)"	"R(Å)" "κ(N/m)" "ReΣ(ps⁻¹)" "ImΣ(ps⁻¹)" "μ(cm²V⁻¹s⁻¹)"]
data1[2:end, 1] = 1.0:400.0
data1[2:end, 2] = ustrip.(rubrene_single_holstein_temp_1to400K.v)
data1[2:end, 3] = ustrip.(rubrene_single_holstein_temp_1to400K.w)
data1[2:end, 4] = ustrip.(rubrene_single_holstein_temp_1to400K.E .|> u"meV")
data1[2:end, 5] = ustrip.(rubrene_mass(rubrene_single_holstein_temp_1to400K))
data1[2:end, 6] = ustrip.(rubrene_radius(rubrene_single_holstein_temp_1to400K))
data1[2:end, 7] = ustrip.(rubrene_spring(rubrene_single_holstein_temp_1to400K))
data1[2:end, 8] = real.(ustrip.(rubrene_memory(rubrene_single_holstein_temp_1to400K)))
data1[2:end, 9] = abs.(imag.(ustrip.(rubrene_memory(rubrene_single_holstein_temp_1to400K))))
data1[2:end, 10] = ustrip.(rubrene_mobility(rubrene_single_holstein_temp_1to400K))
writedlm("data/holstein/variational/rubrene/holstein-rubrene-single-temp-1to400K.dat", data1)
data1
end
  ╠═╡ =#

# ╔═╡ 7552e315-ddf3-4045-a475-3e6b13430d18
md"""
## Multiple Modes
"""

# ╔═╡ da40a936-d2b3-495d-9896-a9775aa165d0
#=╠═╡
begin
data2 = Matrix{Any}(undef, length(rubrene_multi_holstein_temp_1to400K.β) + 1, 10)
data2[1, :] = ["T(K)" "v(THz2π)" "w(THz2π)" "E(meV)" "M(mₑ)"	"R(Å)" "κ(N/m)" "ReΣ(ps⁻¹)" "ImΣ(ps⁻¹)" "μ(cm²V⁻¹s⁻¹)"]
data2[2:end, 1] = 1.0:400.0
data2[2:end, 2] = ustrip.(rubrene_multi_holstein_temp_1to400K.v)
data2[2:end, 3] = ustrip.(rubrene_multi_holstein_temp_1to400K.w)
data2[2:end, 4] = ustrip.(rubrene_multi_holstein_temp_1to400K.E .|> u"meV")
data2[2:end, 5] = ustrip.(rubrene_mass(rubrene_multi_holstein_temp_1to400K))
data2[2:end, 6] = ustrip.(rubrene_radius(rubrene_multi_holstein_temp_1to400K))
data2[2:end, 7] = ustrip.(rubrene_spring(rubrene_multi_holstein_temp_1to400K))
data2[2:end, 8] = real.(ustrip.(rubrene_memory(rubrene_multi_holstein_temp_1to400K)))
data2[2:end, 9] = abs.(imag.(ustrip.(rubrene_memory(rubrene_multi_holstein_temp_1to400K))))
data2[2:end, 10] = ustrip.(rubrene_mobility(rubrene_multi_holstein_temp_1to400K))
writedlm("data/holstein/variational/rubrene/holstein-rubrene-multi-temp-1to400K.dat", data2)
data2
end
  ╠═╡ =#

# ╔═╡ 42d079df-9d76-42f1-bbf1-d5d891f4b048
md"""
# Rubrene | T = 0.0:400.0 K, Ω = 0.01:0.01:30.0 × $rubrene_phonon_frequency
"""

# ╔═╡ 243c71df-b4b3-4ce0-9901-e57b7857dbdc
md"""
## Single Mode
"""

# ╔═╡ 5a35830a-cf26-4426-9201-99667fe2600a
#=╠═╡
begin
data3 = Matrix{Any}(undef, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω) + 1, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.β) + 1)
data3[1, 1] = "ℜΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data3[1, 2:end] = 0.0:400.0
data3[2:end, 1] = ustrip.(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω)
data3[2:end, 2:end] = real.(ustrip.(rubrene_memory(rubrene_single_holstein_temp_0to400K_freq_0to30omega)))'
writedlm("data/holstein/variational/Rubrene/holstein-Rubrene-single-real-memory-temp-0to400K-freq-0to30omega.dat", data3)
data3
end
  ╠═╡ =#

# ╔═╡ 73da4e9d-8556-4bda-ad6f-5b40fe892c27
#=╠═╡
begin
data4 = Matrix{Any}(undef, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω) + 1, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.β) + 1)
data4[1, 1] = "ℑΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data4[1, 2:end] = 0.0:400.0
data4[2:end, 1] = ustrip.(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω)
data4[2:end, 2:end] = abs.(imag.(ustrip.(rubrene_memory(rubrene_single_holstein_temp_0to400K_freq_0to30omega))))'
writedlm("data/holstein/variational/Rubrene/holstein-Rubrene-single-imag-memory-temp-0to400K-freq-0to30omega.dat", data4)
data4
end
  ╠═╡ =#

# ╔═╡ f0384492-ff6f-4767-ba09-ac8b5854cd1f
#=╠═╡
begin
data5 = Matrix{Any}(undef, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω) + 1, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.β) + 1)
data5[1, 1] = "ℜσ(μS)|Ω(THz2π)↓|T(K)→"
data5[1, 2:end] = 0.0:400.0
data5[2:end, 1] = ustrip.(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω)
data5[2:end, 2:end] = ustrip(abs.(real.(rubrene_conductivity(rubrene_single_holstein_temp_0to400K_freq_0to30omega))))
writedlm("data/holstein/variational/Rubrene/holstein-Rubrene-single-real-conductivity-temp-0to400K-freq-0to30omega.dat", data5)
data5
end
  ╠═╡ =#

# ╔═╡ c5f0305e-8601-4cf6-b1d8-9f967ded332d
#=╠═╡
begin
data6 = Matrix{Any}(undef, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω) + 1, length(rubrene_single_holstein_temp_0to400K_freq_0to30omega.β) + 1)
data6[1, 1] = "ℑσ(μS)|Ω(THz2π)↓|T(K)→"
data6[1, 2:end] = 0.0:400.0
data6[2:end, 1] = ustrip.(rubrene_single_holstein_temp_0to400K_freq_0to30omega.Ω)
data6[2:end, 2:end] = ustrip(imag.(rubrene_conductivity(rubrene_single_holstein_temp_0to400K_freq_0to30omega)))
writedlm("data/holstein/variational/Rubrene/holstein-Rubrene-single-imag-conductivity-temp-0to400K-freq-0to30omega.dat", data6)
data6
end
  ╠═╡ =#

# ╔═╡ 84d1260c-f2f5-4ff7-8029-45af59c669e8
md"""
## Multiple Modes
"""

# ╔═╡ 096529e8-2ad0-4cfd-bbfa-de6007d93223
# ╠═╡ disabled = true
#=╠═╡
begin
data7 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data7[1, 1] = "ℜΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data7[1, 2:end] = 0.0:400.0
data7[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data7[2:end, 2:end] = real.(ustrip.(MAPI_memory(MAPIm_frohlich_temp_0to400K_freq_0to30omega)))'
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-real-memory-temp-0to400K-freq-0to30omega.dat", data7)
data7
end
  ╠═╡ =#

# ╔═╡ 01d39122-ceda-4458-9c52-6135fe07afaf
# ╠═╡ disabled = true
#=╠═╡
begin
data8 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data8[1, 1] = "ℑΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data8[1, 2:end] = 0.0:400.0
data8[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data8[2:end, 2:end] = abs.(imag.(ustrip.(MAPI_memory(MAPIm_frohlich_temp_0to400K_freq_0to30omega))))'
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-memory-temp-0to400K-freq-0to30omega.dat", data8)
data8
end
  ╠═╡ =#

# ╔═╡ 05275586-ffc3-4f15-a183-e60e3dcd5beb
# ╠═╡ disabled = true
#=╠═╡
begin
data9 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data9[1, 1] = "ℜσ(μS)|Ω(THz2π)↓|T(K)→"
data9[1, 2:end] = 0.0:400.0
data9[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data9[2:end, 2:end] = ustrip(abs.(real.(MAPI_conductivity(MAPIm_frohlich_temp_0to400K_freq_0to30omega))))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-real-conductivity-temp-0to400K-freq-0to30omega.dat", data9)
data9
end
  ╠═╡ =#

# ╔═╡ 9cfb1ff9-d1fa-4260-b63d-02641fa98fad
# ╠═╡ disabled = true
#=╠═╡
begin
data10 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data10[1, 1] = "ℑσ(μS)|Ω(THz2π)↓|T(K)→"
data10[1, 2:end] = 0.0:400.0
data10[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data10[2:end, 2:end] = ustrip(imag.(MAPI_conductivity(MAPIm_frohlich_temp_0to400K_freq_0to30omega)))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-conductivity-temp-0to400K-freq-0to30omega.dat", data10)
data10
end
  ╠═╡ =#

# ╔═╡ a7826662-223f-46ae-9f2c-8521a6cf220c
#=╠═╡
# 16 threads ~ 8713s
begin
rubrene_multi_holstein_temp_0to400K_freq_0to30omega = holstein(rubrene_multi, LinRange(0,400,401) .* 1u"K", LinRange(0.01,30.0, 3000) .* 
rubrene_phonon_frequency)
save_holstein(rubrene_multi_holstein_temp_0to400K_freq_0to30omega, "data/holstein/variational/Rubrene/holstein-rubrene-multi-temp-0to400K-freq-0to30omega")
end
  ╠═╡ =#

# ╔═╡ 0c48545a-1646-4405-9e9f-6e53d04cc0a7
# ╠═╡ disabled = true
#=╠═╡
rubrene_multi_holstein_temp_0to400K_freq_0to30omega = load_holstein("data/holstein/variational/Rubrene/holstein-rubrene-multi-temp-0to400K-freq-0to30omega.jld")
  ╠═╡ =#

# ╔═╡ c670dd62-a06d-4fb4-83d7-ae5fdb54baf3
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 8713s
begin
rubrene_single_holstein_temp_0to400K_freq_0to30omega = holstein(rubrene_single, LinRange(0,400,401) .* 1u"K", LinRange(0.01,30.0, 3000) .* 
rubrene_phonon_frequency)
save_holstein(rubrene_single_holstein_temp_0to400K_freq_0to30omega, "data/holstein/variational/Rubrene/holstein-rubrene-single-temp-0to400K-freq-0to30omega")
end
  ╠═╡ =#

# ╔═╡ 76efc24c-812a-4e22-86e1-0f928ef13677
# ╠═╡ disabled = true
#=╠═╡
begin
rubrene_single = HolsteinMaterial("C18H8[C6H5]4", rubrene_transfer_integral, rubrene_phonon_frequency, sqrt.(rubrene_holstein_coupling_energy .^2 .+ rubrene_peierls_coupling_energy .^2), rubrene_lattice_constant, 1.0)
save_holstein_material(rubrene_single, "data/holstein/materials/Rubrene/rubrene-single")
end
  ╠═╡ =#

# ╔═╡ 39ad223a-76e2-46f6-9ff9-b1e206b8bca1
#=╠═╡
rubrene_single_holstein_temp_1to400K = load_holstein("data/holstein/variational/Rubrene/holstein-rubrene-single-temp-1to400K.jld")
  ╠═╡ =#

# ╔═╡ 0c1b1590-6fed-4cf5-b829-c65fdfccf166
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 490s
begin
rubrene_multi_holstein_temp_1to400K = holstein(rubrene_multi, LinRange(1,400,400) .* 1u"K", eps())
save_holstein(rubrene_multi_holstein_temp_1to400K, "data/holstein/variational/Rubrene/holstein-rubrene-multi-temp-1to400K")
end
  ╠═╡ =#

# ╔═╡ 60bdbb74-e4d5-4d6e-b374-45b0ead3dcf9
#=╠═╡
rubrene_single_holstein_temp_0to400K_freq_0to30omega = load_holstein("data/holstein/variational/Rubrene/holstein-rubrene-single-temp-0to400K-freq-0to30omega.jld")
  ╠═╡ =#

# ╔═╡ cff673c1-b106-44b7-89ed-4d3a3cd588a3
# ╠═╡ disabled = true
#=╠═╡
begin
rubrene_multi = HolsteinMaterial("C18H8[C6H5]4", rubrene_transfer_integral, rubrene_phonon_frequencies, sqrt.(rubrene_holstein_coupling_energies .^2 .+ rubrene_peierls_coupling_energies .^2), rubrene_lattice_constant, 1.0)
save_holstein_material(rubrene_multi, "data/holstein/materials/Rubrene/rubrene-multi")
end
  ╠═╡ =#

# ╔═╡ e39f0dca-d015-42f1-a49d-ee88b938e545
#=╠═╡
rubrene_multi = load_holstein_material("data/holstein/materials/Rubrene/rubrene-multi.jld")
  ╠═╡ =#

# ╔═╡ 568d6894-2895-4c43-985e-32e0fad66210
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 273s
begin
rubrene_single_holstein_temp_1to400K = holstein(rubrene_single, LinRange(1,400,400) .* 1u"K", eps())
save_holstein(rubrene_single_holstein_temp_1to400K, "data/holstein/variational/Rubrene/holstein-rubrene-single-temp-1to400K")
end
  ╠═╡ =#

# ╔═╡ 8ddf0c8d-794f-4188-b38f-0dea63dd4656
#=╠═╡
rubrene_multi_holstein_temp_1to400K = load_holstein("data/holstein/variational/Rubrene/holstein-rubrene-multi-temp-1to400K.jld")
  ╠═╡ =#

# ╔═╡ 9bcc691e-e704-464a-85a5-fd6bcdbecc6d
#=╠═╡
rubrene_single = load_holstein_material("data/holstein/materials/Rubrene/rubrene-single.jld")
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─69067b24-4916-11ef-037f-0902e798e948
# ╟─ad02f2ba-b5d6-44ba-afcd-279be120007a
# ╟─50109bde-f4ee-4918-9c74-753cad69a2e6
# ╟─bca41c3b-7daa-4445-997e-b3fef2e67d56
# ╟─e6f4c89c-32cb-419f-9a92-fd4837632b81
# ╟─6d58793a-126c-44d5-adc7-b4acfbdae12a
# ╟─cbf984e6-602b-4753-9ecc-db2ddbe28a6b
# ╟─0135db00-7d9a-45ec-b9b9-bd4359a02089
# ╟─1f1e81eb-a571-46c2-8618-e252cf69d205
# ╟─799e2865-eb50-49ca-8aaf-e1d648ae2b25
# ╟─17068899-f58f-4415-9057-06531b35ab66
# ╟─762f3da9-af65-4b10-8210-6c4255ad9526
# ╟─82184ac9-eebd-4a23-8347-965206153a16
# ╟─fe1b43ac-be1c-4a77-bee3-dc3d5b227e2e
# ╟─c46076de-52ee-48ed-ab68-e07eefc168f3
# ╟─662636a3-c9b4-4552-96ce-360d3e0fc93d
# ╟─2c9625e7-11e5-4d56-84e8-3d15ba0a6c6a
# ╟─0ba581f6-1e4c-4960-8668-c00c1084bd64
# ╟─5f8a7623-bcb4-4797-98a9-d26450d39409
# ╟─6298f629-457d-4881-82c4-2c356de25079
# ╟─78b458f9-8c59-4949-9555-64a1eda46709
# ╟─9060c1ad-5ce1-4bd5-a2af-45af9c04f0c3
# ╟─d71db010-e118-4f22-99a3-8327e1b0e303
# ╟─a1d005ec-94cf-48f4-a667-ea9733497bab
# ╟─3851a3fe-6f3b-4099-b7a3-038dd1f5bf63
# ╟─e74f64c8-bedc-461e-9cf0-407f8f1d4313
# ╟─76efc24c-812a-4e22-86e1-0f928ef13677
# ╟─9bcc691e-e704-464a-85a5-fd6bcdbecc6d
# ╟─39eda869-dfab-47e7-9c8d-2ca734eed33e
# ╟─cff673c1-b106-44b7-89ed-4d3a3cd588a3
# ╟─e39f0dca-d015-42f1-a49d-ee88b938e545
# ╟─5270ea61-c216-47e8-a9de-87635dabecec
# ╟─c200c995-4eb8-4347-82fc-7dc170401df8
# ╟─568d6894-2895-4c43-985e-32e0fad66210
# ╟─39ad223a-76e2-46f6-9ff9-b1e206b8bca1
# ╟─8acd1416-755e-405f-852f-884a7976ca82
# ╟─7552e315-ddf3-4045-a475-3e6b13430d18
# ╟─0c1b1590-6fed-4cf5-b829-c65fdfccf166
# ╟─8ddf0c8d-794f-4188-b38f-0dea63dd4656
# ╟─da40a936-d2b3-495d-9896-a9775aa165d0
# ╟─42d079df-9d76-42f1-bbf1-d5d891f4b048
# ╟─243c71df-b4b3-4ce0-9901-e57b7857dbdc
# ╟─c670dd62-a06d-4fb4-83d7-ae5fdb54baf3
# ╟─60bdbb74-e4d5-4d6e-b374-45b0ead3dcf9
# ╟─5a35830a-cf26-4426-9201-99667fe2600a
# ╟─73da4e9d-8556-4bda-ad6f-5b40fe892c27
# ╟─f0384492-ff6f-4767-ba09-ac8b5854cd1f
# ╟─c5f0305e-8601-4cf6-b1d8-9f967ded332d
# ╟─84d1260c-f2f5-4ff7-8029-45af59c669e8
# ╠═a7826662-223f-46ae-9f2c-8521a6cf220c
# ╠═0c48545a-1646-4405-9e9f-6e53d04cc0a7
# ╟─096529e8-2ad0-4cfd-bbfa-de6007d93223
# ╟─01d39122-ceda-4458-9c52-6135fe07afaf
# ╟─05275586-ffc3-4f15-a183-e60e3dcd5beb
# ╟─9cfb1ff9-d1fa-4260-b63d-02641fa98fad
