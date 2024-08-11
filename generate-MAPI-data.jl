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
MAPI_memory(frohlich) = frohlich.Σ .|> u"ps^-1"

# ╔═╡ 17068899-f58f-4415-9057-06531b35ab66
md"""
# MAPI Literature Data
"""

# ╔═╡ f7baaa29-046b-41c6-aa95-18a4591555ff
MAPI = [
   # 96.20813558773261 0.4996300522819191
   # 93.13630357703363 1.7139631746083817
   # 92.87834578121567 0.60108592692181
   # 92.4847918585963 0.0058228799414729
   # 92.26701437594754 0.100590086574602
   # 89.43972834606603 0.006278895133832249
   # 46.89209141511332 0.2460894564364346
   # 46.420949316788 0.14174282581124137
   # 44.0380222871706 0.1987196948553428
   # 42.89702947649343 0.011159939465770681
   # 42.67180170168193 0.02557751102757614
   # 41.46971205834201 0.012555230726601503
   # 37.08982543385215 0.00107488277468418
   # 36.53555265689563 0.02126940080871224
   # 30.20608114002676 0.009019481779712388
   # 27.374810898415028 0.03994453721421388
   # 26.363055017011728 0.05011922682554448
   # 9.522966890022039 0.00075631870522737
   4.016471586720514 0.08168931020200264
   3.887605410774121 0.006311654262282101
   3.5313112232401513 0.05353548710183397
   2.755392921480459 0.021303020776321225
   2.4380741812443247 0.23162784335484837
   2.2490917637719408 0.2622203718355982
   2.079632190634424 0.23382298607799906
   2.0336707697261187 0.0623239656843172
   1.5673011873879714 0.0367465760261409
   1.0188379384951798 0.0126328938653956
   1.0022960504442775 0.006817361620021601
   0.9970130778462072 0.0103757951973341
   0.9201781906386209 0.01095811116040592
   0.800604081794174 0.0016830270365341532
   0.5738689505255512 0.00646428491253749
   # 0.022939578929507105 8.355742795827834e-05   # Acoustic modes!
   # 0.04882611767873102 8.309858592685e-06
   # 0.07575149723846182 2.778248540373041e-05
]

# ╔═╡ 6909665a-1496-4b01-88eb-dbf499c953ae
begin
	phonon_frequencies = MAPI[:, 1] .* u"THz2π"
	phonon_frequency = 1.62289245u"THz2π"
	band_mass = 0.12u"me"
	infrared_activities = MAPI[:, 2] .* u"q^2 / u"
	unitcell_volume = (6.29u"Å")^3
	electronic_dielectric = 4.5u"ϵ0"
	ionic_dielectric = ϵ_ionic_mode.(phonon_frequencies, infrared_activities, unitcell_volume)
	static_dielectric = 4.5u"ϵ0" + sum(ionic_dielectric)
end

# ╔═╡ 3851a3fe-6f3b-4099-b7a3-038dd1f5bf63
md"""
# MAPI Material Data
"""

# ╔═╡ e74f64c8-bedc-461e-9cf0-407f8f1d4313
md"""
## Single Mode
"""

# ╔═╡ 76efc24c-812a-4e22-86e1-0f928ef13677
# ╠═╡ disabled = true
#=╠═╡
begin
MAPIs = frohlichmaterial("[CH3NH3]PbX3", static_dielectric, electronic_dielectric, phonon_frequency, band_mass, unitcell_volume)
save_frohlich_material(MAPIs, "data/frohlich/materials/MAPI/MAPI-single")
end
  ╠═╡ =#

# ╔═╡ 9bcc691e-e704-464a-85a5-fd6bcdbecc6d
MAPIs = load_frohlich_material("data/frohlich/materials/MAPI/MAPI-single.jld")

# ╔═╡ e6f4c89c-32cb-419f-9a92-fd4837632b81
MAPI_mass(frohlich) = frohlich.v.^2 ./ frohlich.w.^2 .* MAPIs.mb

# ╔═╡ 6d58793a-126c-44d5-adc7-b4acfbdae12a
MAPI_spring(frohlich) = (frohlich.v.^2 .- frohlich.w.^2) .* MAPIs.mb .|> u"N/m"

# ╔═╡ 78dc4361-a01e-4aed-880f-65d3f7e75e9e
MAPI_radius(frohlich; dims=3) = dims / 2 .* frohlich.v ./ (frohlich.v.^2 .- frohlich.w.^2) .* coth.(pustrip.(frohlich.β .* frohlich.v ./ 2)) .* sqrt(ħ_pu * MAPIs.ω_LO / MAPIs.mb) .|> u"Å"

# ╔═╡ 7eae6361-48bf-4994-a529-0cb4b9c84e0b
MAPI_radius2(frohlich; dims=3) = dims / 2 .* pustrip.(frohlich.v[:,1,:])' ./ (pustrip.(frohlich.v[:,1,:])'.^2 .- pustrip.(frohlich.w[:,1,:])'.^2) .* coth.(pustrip.(frohlich.β .* pustrip.(frohlich.v[:,1,:])' ./ 2)) .* sqrt(ħ_pu / MAPIs.mb / MAPIs.ω_LO) .|> u"Å"

# ╔═╡ 0135db00-7d9a-45ec-b9b9-bd4359a02089
MAPI_mobility(frohlich) = abs.(1 ./ imag.(frohlich.Σ)) * e_pu / MAPIs.mb .|> u"cm^2/V/s"

# ╔═╡ 799e2865-eb50-49ca-8aaf-e1d648ae2b25
MAPI_conductivity(frohlich) = -im ./ (frohlich.Ω .+ frohlich.Σ') .* e_pu^2 / MAPIs.mb / (ħ_pu / MAPIs.mb / MAPIs.ω_LO) .|> u"μS"

# ╔═╡ 39eda869-dfab-47e7-9c8d-2ca734eed33e
md"""
## Multiple Modes
"""

# ╔═╡ cff673c1-b106-44b7-89ed-4d3a3cd588a3
# ╠═╡ disabled = true
#=╠═╡
begin
MAPIm = frohlichmaterial("[CH3NH3]PbX3", static_dielectric, electronic_dielectric, ionic_dielectric, phonon_frequencies, band_mass, unitcell_volume)
save_frohlich_material(MAPIm, "data/frohlich/materials/MAPI/MAPI-multi")
end
  ╠═╡ =#

# ╔═╡ e39f0dca-d015-42f1-a49d-ee88b938e545
MAPIm = load_frohlich_material("data/frohlich/materials/MAPI/MAPI-multi.jld")

# ╔═╡ 5270ea61-c216-47e8-a9de-87635dabecec
md"""
# MAPI | T = 1:400 K 
"""

# ╔═╡ c200c995-4eb8-4347-82fc-7dc170401df8
md"""
## Single Mode
"""

# ╔═╡ 568d6894-2895-4c43-985e-32e0fad66210
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 27s
begin
MAPIs_frohlich_temp_1to400K = frohlich(MAPIs, LinRange(1,400,400) .* 1u"K", eps())
save_frohlich(MAPIs_frohlich_temp_1to400K, "data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K")
end
  ╠═╡ =#

# ╔═╡ 39ad223a-76e2-46f6-9ff9-b1e206b8bca1
MAPIs_frohlich_temp_1to400K = load_frohlich("data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.jld")

# ╔═╡ afd07d59-b067-4350-80bf-e26cc3e2c086
begin
data1 = Matrix{Any}(undef, length(MAPIs_frohlich_temp_1to400K.β) + 1, 10)
data1[1, :] = ["T(K)" "v(THz2π)" "w(THz2π)" "E(meV)" "M(mₑ)"	"R(Å)" "κ(N/m)" "ReΣ(ps⁻¹)" "ImΣ(ps⁻¹)" "μ(cm²V⁻¹s⁻¹)"]
data1[2:end, 1] = 1.0:400.0
data1[2:end, 2] = ustrip.(MAPIs_frohlich_temp_1to400K.v)
data1[2:end, 3] = ustrip.(MAPIs_frohlich_temp_1to400K.w)
data1[2:end, 4] = ustrip.(MAPIs_frohlich_temp_1to400K.E .|> u"meV")
data1[2:end, 5] = ustrip.(MAPI_mass(MAPIs_frohlich_temp_1to400K))
data1[2:end, 6] = ustrip.(MAPI_radius2(MAPIs_frohlich_temp_1to400K))
data1[2:end, 7] = ustrip.(MAPI_spring(MAPIs_frohlich_temp_1to400K))
data1[2:end, 8] = real.(ustrip.(MAPI_memory(MAPIs_frohlich_temp_1to400K)))
data1[2:end, 9] = abs.(imag.(ustrip.(MAPI_memory(MAPIs_frohlich_temp_1to400K))))
data1[2:end, 10] = ustrip.(MAPI_mobility(MAPIs_frohlich_temp_1to400K))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.dat", data1)
data1
end

# ╔═╡ 7552e315-ddf3-4045-a475-3e6b13430d18
md"""
## Multiple Modes
"""

# ╔═╡ 0c1b1590-6fed-4cf5-b829-c65fdfccf166
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 490s
begin
MAPIm_frohlich_temp_1to400K = frohlich(MAPIm, LinRange(1,400,400) .* 1u"K", eps())
save_frohlich(MAPIm_frohlich_temp_1to400K, "data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K")
end
  ╠═╡ =#

# ╔═╡ 8ddf0c8d-794f-4188-b38f-0dea63dd4656
MAPIm_frohlich_temp_1to400K = load_frohlich("data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.jld")

# ╔═╡ da40a936-d2b3-495d-9896-a9775aa165d0
begin
data2 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_1to400K.β) + 1, 10)
data2[1, :] = ["T(K)" "v(THz2π)" "w(THz2π)" "E(meV)" "M(mₑ)"	"R(Å)" "κ(N/m)" "ReΣ(ps⁻¹)" "ImΣ(ps⁻¹)" "μ(cm²V⁻¹s⁻¹)"]
data2[2:end, 1] = 1.0:400.0
data2[2:end, 2] = ustrip.(MAPIm_frohlich_temp_1to400K.v)
data2[2:end, 3] = ustrip.(MAPIm_frohlich_temp_1to400K.w)
data2[2:end, 4] = ustrip.(MAPIm_frohlich_temp_1to400K.E .|> u"meV")
data2[2:end, 5] = ustrip.(MAPI_mass(MAPIm_frohlich_temp_1to400K))
data2[2:end, 6] = ustrip.(MAPI_radius(MAPIm_frohlich_temp_1to400K))
data2[2:end, 7] = ustrip.(MAPI_spring(MAPIm_frohlich_temp_1to400K))
data2[2:end, 8] = real.(ustrip.(MAPI_memory(MAPIm_frohlich_temp_1to400K)))
data2[2:end, 9] = abs.(imag.(ustrip.(MAPI_memory(MAPIm_frohlich_temp_1to400K))))
data2[2:end, 10] = ustrip.(MAPI_mobility(MAPIm_frohlich_temp_1to400K))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.dat", data2)
data2
end

# ╔═╡ 42d079df-9d76-42f1-bbf1-d5d891f4b048
md"""
# MAPI | T = 0.0:400.0 K, Ω = 0.01:0.01:30.0 × $phonon_frequency
"""

# ╔═╡ 243c71df-b4b3-4ce0-9901-e57b7857dbdc
md"""
## Single Mode
"""

# ╔═╡ c670dd62-a06d-4fb4-83d7-ae5fdb54baf3
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 13520s
begin
MAPIs_frohlich_temp_0to400K_freq_0to30omega = frohlich(MAPIs, LinRange(0,400,401) .* 1u"K", LinRange(0.01,30.0, 3000) .* phonon_frequency)
save_frohlich(MAPIs_frohlich_temp_0to400K_freq_0to30omega, "data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-0to400K-freq-0to30omega")
end
  ╠═╡ =#

# ╔═╡ 60bdbb74-e4d5-4d6e-b374-45b0ead3dcf9
MAPIs_frohlich_temp_0to400K_freq_0to30omega = load_frohlich("data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-0to400K-freq-0to30omega.jld")

# ╔═╡ 5a35830a-cf26-4426-9201-99667fe2600a
begin
data3 = Matrix{Any}(undef, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data3[1, 1] = "ℜΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data3[1, 2:end] = 0.0:400.0
data3[2:end, 1] = ustrip.(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω)
data3[2:end, 2:end] = real.(ustrip.(MAPI_memory(MAPIs_frohlich_temp_0to400K_freq_0to30omega)))'
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-single-real-memory-temp-0to400K-freq-0to30omega.dat", data3)
data3
end

# ╔═╡ 73da4e9d-8556-4bda-ad6f-5b40fe892c27
begin
data4 = Matrix{Any}(undef, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data4[1, 1] = "ℑΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data4[1, 2:end] = 0.0:400.0
data4[2:end, 1] = ustrip.(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω)
data4[2:end, 2:end] = abs.(imag.(ustrip.(MAPI_memory(MAPIs_frohlich_temp_0to400K_freq_0to30omega))))'
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-single-imag-memory-temp-0to400K-freq-0to30omega.dat", data4)
data4
end

# ╔═╡ f0384492-ff6f-4767-ba09-ac8b5854cd1f
begin
data5 = Matrix{Any}(undef, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data5[1, 1] = "ℜσ(μS)|Ω(THz2π)↓|T(K)→"
data5[1, 2:end] = 0.0:400.0
data5[2:end, 1] = ustrip.(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω)
data5[2:end, 2:end] = ustrip(abs.(real.(MAPI_conductivity(MAPIs_frohlich_temp_0to400K_freq_0to30omega))))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0to400K-freq-0to30omega.dat", data5)
data5
end

# ╔═╡ c5f0305e-8601-4cf6-b1d8-9f967ded332d
begin
data6 = Matrix{Any}(undef, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIs_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data6[1, 1] = "ℑσ(μS)|Ω(THz2π)↓|T(K)→"
data6[1, 2:end] = 0.0:400.0
data6[2:end, 1] = ustrip.(MAPIs_frohlich_temp_0to400K_freq_0to30omega.Ω)
data6[2:end, 2:end] = ustrip(imag.(MAPI_conductivity(MAPIs_frohlich_temp_0to400K_freq_0to30omega)))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-single-imag-conductivity-temp-0to400K-freq-0to30omega.dat", data6)
data6
end

# ╔═╡ 84d1260c-f2f5-4ff7-8029-45af59c669e8
md"""
## Multiple Modes
"""

# ╔═╡ a7826662-223f-46ae-9f2c-8521a6cf220c
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 26348s
begin
MAPIm_frohlich_temp_0to400K_freq_0to30omega = frohlich(MAPIm, LinRange(0,400,401) .* 1u"K", LinRange(0.01,30.0, 3000) .* phonon_frequency)
save_frohlich(MAPIm_frohlich_temp_0to400K_freq_0to30omega, "data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-0to400K-freq-0to30omega")
end
  ╠═╡ =#

# ╔═╡ 0c48545a-1646-4405-9e9f-6e53d04cc0a7
MAPIm_frohlich_temp_0to400K_freq_0to30omega = load_frohlich("data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-0to400K-freq-0to30omega.jld")

# ╔═╡ 096529e8-2ad0-4cfd-bbfa-de6007d93223
begin
data7 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data7[1, 1] = "ℜΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data7[1, 2:end] = 0.0:400.0
data7[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data7[2:end, 2:end] = real.(ustrip.(MAPI_memory(MAPIm_frohlich_temp_0to400K_freq_0to30omega)))'
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-real-memory-temp-0to400K-freq-0to30omega.dat", data7)
data7
end

# ╔═╡ 01d39122-ceda-4458-9c52-6135fe07afaf
begin
data8 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data8[1, 1] = "ℑΣ(ps⁻¹)|Ω(THz2π)↓|T(K)→"
data8[1, 2:end] = 0.0:400.0
data8[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data8[2:end, 2:end] = abs.(imag.(ustrip.(MAPI_memory(MAPIm_frohlich_temp_0to400K_freq_0to30omega))))'
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-memory-temp-0to400K-freq-0to30omega.dat", data8)
data8
end

# ╔═╡ 05275586-ffc3-4f15-a183-e60e3dcd5beb
begin
data9 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data9[1, 1] = "ℜσ(μS)|Ω(THz2π)↓|T(K)→"
data9[1, 2:end] = 0.0:400.0
data9[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data9[2:end, 2:end] = ustrip(abs.(real.(MAPI_conductivity(MAPIm_frohlich_temp_0to400K_freq_0to30omega))))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-real-conductivity-temp-0to400K-freq-0to30omega.dat", data9)
data9
end

# ╔═╡ 9cfb1ff9-d1fa-4260-b63d-02641fa98fad
begin
data10 = Matrix{Any}(undef, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω) + 1, length(MAPIm_frohlich_temp_0to400K_freq_0to30omega.β) + 1)
data10[1, 1] = "ℑσ(μS)|Ω(THz2π)↓|T(K)→"
data10[1, 2:end] = 0.0:400.0
data10[2:end, 1] = ustrip.(MAPIm_frohlich_temp_0to400K_freq_0to30omega.Ω)
data10[2:end, 2:end] = ustrip(imag.(MAPI_conductivity(MAPIm_frohlich_temp_0to400K_freq_0to30omega)))
writedlm("data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-conductivity-temp-0to400K-freq-0to30omega.dat", data10)
data10
end

# ╔═╡ Cell order:
# ╟─69067b24-4916-11ef-037f-0902e798e948
# ╟─ad02f2ba-b5d6-44ba-afcd-279be120007a
# ╟─50109bde-f4ee-4918-9c74-753cad69a2e6
# ╟─bca41c3b-7daa-4445-997e-b3fef2e67d56
# ╟─e6f4c89c-32cb-419f-9a92-fd4837632b81
# ╟─6d58793a-126c-44d5-adc7-b4acfbdae12a
# ╟─78dc4361-a01e-4aed-880f-65d3f7e75e9e
# ╟─7eae6361-48bf-4994-a529-0cb4b9c84e0b
# ╟─0135db00-7d9a-45ec-b9b9-bd4359a02089
# ╟─1f1e81eb-a571-46c2-8618-e252cf69d205
# ╟─799e2865-eb50-49ca-8aaf-e1d648ae2b25
# ╟─17068899-f58f-4415-9057-06531b35ab66
# ╟─f7baaa29-046b-41c6-aa95-18a4591555ff
# ╠═6909665a-1496-4b01-88eb-dbf499c953ae
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
# ╟─afd07d59-b067-4350-80bf-e26cc3e2c086
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
# ╟─a7826662-223f-46ae-9f2c-8521a6cf220c
# ╟─0c48545a-1646-4405-9e9f-6e53d04cc0a7
# ╟─096529e8-2ad0-4cfd-bbfa-de6007d93223
# ╟─01d39122-ceda-4458-9c52-6135fe07afaf
# ╟─05275586-ffc3-4f15-a183-e60e3dcd5beb
# ╟─9cfb1ff9-d1fa-4260-b63d-02641fa98fad
