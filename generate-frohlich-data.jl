### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bfd6e85c-48e2-11ef-2751-3147edce85d8
using Pkg, Unitful, Plots, DelimitedFiles

# ╔═╡ 68c1b156-ee71-4717-94d6-914500fab1ff
Pkg.activate(".")

# ╔═╡ 52e5925c-a82a-4799-97b4-a96a0798fa9b
using PolaronThesis

# ╔═╡ b12c5fa7-a588-4ce1-a4b5-1cf598380903
md"""
# Polaron Equations
"""

# ╔═╡ 4d74ec93-b41d-40af-b8e6-0652b2f10faf
polaron_mass(frohlich) = frohlich.v.^2 ./ frohlich.w.^2

# ╔═╡ 83427388-3d3f-4f05-867a-6984454cfdf9
polaron_spring(frohlich) = frohlich.v.^2 .- frohlich.w.^2

# ╔═╡ 8564e569-76ed-439a-aa13-291538cfff13
polaron_radius(frohlich; dims=3) = dims / 2 .* frohlich.v ./ (frohlich.v.^2 .- frohlich.w.^2) .* coth.(pustrip.(frohlich.β .* frohlich.v ./ 2))

# ╔═╡ 828f01c9-5726-4c8a-a5b1-bd6e83f29662
polaron_radius2(frohlich; dims=3) = dims / 2 .* pustrip.(frohlich.v[:,1,:])' ./ (pustrip.(frohlich.v[:,1,:])'.^2 .- pustrip.(frohlich.w[:,1,:])'.^2) .* coth.(pustrip.(frohlich.β .* pustrip.(frohlich.v[:,1,:])' ./ 2))

# ╔═╡ 6088c544-f44b-49c1-a503-ca881418e2f9
polaron_mobility(frohlich) = abs.(1 ./ imag.(frohlich.Σ))

# ╔═╡ 8420455b-5772-413d-8e19-4e541e379142
polaron_conductivity(frohlich) = -im ./ (frohlich.Ω .+ frohlich.Σ')

# ╔═╡ e03e239c-a77d-4591-aaa0-e02ce4cc0292
md"""
# 2D and 3D Frohlich α = 0.1:0.1:12: v ω₀, w ω₀, E ħω₀, M m₀, R √(ħ/m₀ω₀), κ m₀ω₀²
"""

# ╔═╡ 9d95d5b9-1eee-4d7c-acd2-a2177f620fb7
md"""
## 2D
"""

# ╔═╡ 174468b5-1853-4eeb-858f-3a65f9c8a822
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 1s
begin
frohlich_2d_alpha_0to12_beta_inf = frohlich(LinRange(0.1,12,120), 1.0; dims = 2)
save_frohlich(frohlich_2d_alpha_0to12_beta_inf, "data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf")
end
  ╠═╡ =#

# ╔═╡ c19ed360-aa81-4ed1-a5b6-0359b01e9a5f
frohlich_2d_alpha_0to12_beta_inf = load_frohlich("data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.jld")

# ╔═╡ 2cb035de-37aa-434b-b61e-a9106bac6b86
begin
data1 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_inf.α) + 1, 7)
data1[1, :] = ["α" "v(ω₀)" "w(ω₀)" "E(ħω₀)" "M(m₀)"	"R(r₀)" "κ(m₀ω₀²)"]
data1[2:end, 1] = frohlich_2d_alpha_0to12_beta_inf.α
data1[2:end, 2] = ustrip.(frohlich_2d_alpha_0to12_beta_inf.v)
data1[2:end, 3] = ustrip.(frohlich_2d_alpha_0to12_beta_inf.w)
data1[2:end, 4] = ustrip.(frohlich_2d_alpha_0to12_beta_inf.E)
data1[2:end, 5] = ustrip.(polaron_mass(frohlich_2d_alpha_0to12_beta_inf))
data1[2:end, 6] = ustrip.(polaron_radius(frohlich_2d_alpha_0to12_beta_inf))
data1[2:end, 7] = ustrip.(polaron_spring(frohlich_2d_alpha_0to12_beta_inf))
writedlm("data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.dat", data1)
data1
end

# ╔═╡ 90a6cf36-3b99-4841-bb6a-5c0192c11c37
md"""
## 3D
"""

# ╔═╡ 05caea00-de80-4357-83b3-8d884b8e0573
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 1s
begin
frohlich_3d_alpha_0to12_beta_inf = frohlich(LinRange(0.1,12,120), 1.0)
save_frohlich(frohlich_3d_alpha_0to12_beta_inf, "data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf")
end
  ╠═╡ =#

# ╔═╡ 9a91d18a-b49f-46b5-96e1-f577b142ceac
frohlich_3d_alpha_0to12_beta_inf = load_frohlich("data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.jld")

# ╔═╡ 477eaae2-731a-439a-9888-1d00f42943a8
begin
data2 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_inf.α) + 1, 7)
data2[1, :] = ["α" "v(ω₀)" "w(ω₀)" "E(ħω₀)" "M(m₀)"	"R(r₀)" "κ(m₀ω₀²)"]
data2[2:end, 1] = frohlich_3d_alpha_0to12_beta_inf.α
data2[2:end, 2] = ustrip.(frohlich_3d_alpha_0to12_beta_inf.v)
data2[2:end, 3] = ustrip.(frohlich_3d_alpha_0to12_beta_inf.w)
data2[2:end, 4] = ustrip.(frohlich_3d_alpha_0to12_beta_inf.E)
data2[2:end, 5] = ustrip.(polaron_mass(frohlich_3d_alpha_0to12_beta_inf))
data2[2:end, 6] = ustrip.(polaron_radius(frohlich_3d_alpha_0to12_beta_inf))
data2[2:end, 7] = ustrip.(polaron_spring(frohlich_3d_alpha_0to12_beta_inf))
writedlm("data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat", data2)
data2
end

# ╔═╡ 92b52583-b5db-43cd-8ee9-06098deece54
md"""
# 2D and 3D Frohlich α = 0.1:0.1:12, T = LogRange(0.0625,32,128) ħω₀/kB
"""

# ╔═╡ 4589e7b4-5bfe-437f-927b-35626eb52b84
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 4475s
begin
frohlich_2d_alpha_0to12_beta_003125to16 = frohlich(LinRange(12,0.1,120), 1.0, LogRange(16.0, 0.03125, 127), eps(); dims = 2)
save_frohlich(frohlich_2d_alpha_0to12_beta_003125to16, "data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-003125to16")
end
  ╠═╡ =#

# ╔═╡ ae307e8c-13c4-4d5b-8f99-4c7242df1565
frohlich_2d_alpha_0to12_beta_003125to16 = load_frohlich("data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-003125to16.jld")

# ╔═╡ 808e75f2-7956-4eef-bec8-0564c8750aec
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 4406s
begin
frohlich_3d_alpha_0to12_beta_003125to16 = frohlich(LinRange(12,0.1,120), 1.0, LogRange(16.0, 0.03125, 127), eps())
save_frohlich(frohlich_3d_alpha_0to12_beta_003125to16, "data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-003125to16")
end
  ╠═╡ =#

# ╔═╡ a9de638c-c8f7-438a-89fc-e362d7e725cb
frohlich_3d_alpha_0to12_beta_003125to16 = load_frohlich("data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-003125to16.jld")

# ╔═╡ 8810d4e6-4ce8-4c29-af10-e34d7470a062
md"""
## Energy ħω₀
"""

# ╔═╡ 6a881a63-e422-4389-b25b-23a2c6d82d29
md"""
### 2D
"""

# ╔═╡ c8d1b9fc-9c07-4ebb-8694-05f9586943eb
begin
data3 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data3[1, 1] = "E(ħω₀)|T(ħω₀/kB)↓|α→"
data3[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data3[1, 2:end] = 0.1:0.1:12
data3[2:end, 2:end] = reverse(ustrip(frohlich_2d_alpha_0to12_beta_003125to16.E[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat", data3)
data3
end

# ╔═╡ f134f1d4-33c1-48c8-89ec-0f9f2d457e3e
md"""
### 3D
"""

# ╔═╡ 99faadff-2dd4-46ce-9a8c-d8faf850c9b7
begin
data4 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_3d_alpha_0to12_beta_003125to16.α) + 1)
data4[1, 1] = "E(ħω₀)|T(ħω₀/kB)↓|α→"
data4[2:end, 1] = ustrip(1 ./ frohlich_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data4[1, 2:end] = 0.1:0.1:12
data4[2:end, 2:end] = reverse(ustrip(frohlich_3d_alpha_0to12_beta_003125to16.E[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-E-alpha-0to12-temp-00625to32.dat", data4)
data4
end

# ╔═╡ 958e1299-defb-410d-ac85-c176f621462c
md"""
## v ω₀
"""

# ╔═╡ caa54305-0f13-47a4-850d-a3b14d975ad1
md"""
### 2D
"""

# ╔═╡ e5d9f83f-1c5d-4e0a-9842-554846c227dc
begin
data5 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data5[1, 1] = "v(ω₀)|T(ħω₀/kB)↓|α→"
data5[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data5[1, 2:end] = 0.1:0.1:12
data5[2:end, 2:end] = reverse(ustrip(frohlich_2d_alpha_0to12_beta_003125to16.v[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-v-alpha-0to12-temp-00625to32.dat", data5)
data5
end

# ╔═╡ 13b93723-dcd3-4bc5-9935-30635ec2bbab
md"""
### 3D
"""

# ╔═╡ 1e306f70-be2b-4704-aa70-9d94d85d8859
begin
data6 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_3d_alpha_0to12_beta_003125to16.α) + 1)
data6[1, 1] = "v(ω₀)|T(ħω₀/kB)↓|α→"
data6[2:end, 1] = ustrip(1 ./ frohlich_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data6[1, 2:end] = 0.1:0.1:12
data6[2:end, 2:end] = reverse(ustrip(frohlich_3d_alpha_0to12_beta_003125to16.v[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-v-alpha-0to12-temp-00625to32.dat", data6)
data6
end

# ╔═╡ c417194c-084c-440e-9bed-a5d7b9219ac8
md"""
## w ω₀
"""

# ╔═╡ 2f330d5d-ebd5-4e56-b9d6-5e4dc8c48ac4
md"""
### 2D
"""

# ╔═╡ 670e50ac-afea-4e62-aef7-9bd12a8541d5
begin
data7 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data7[1, 1] = "w(ω₀)|T(ħω₀/kB)↓|α→"
data7[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data7[1, 2:end] = 0.1:0.1:12
data7[2:end, 2:end] = reverse(ustrip(frohlich_2d_alpha_0to12_beta_003125to16.w[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-w-alpha-0to12-temp-00625to32.dat", data7)
data7
end

# ╔═╡ e8c3adce-01b9-4ffe-a946-585cce2e4ae4
md"""
### 3D
"""

# ╔═╡ c3a3a96a-cc7c-4f13-9a01-14473007fcd9
begin
data8 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_3d_alpha_0to12_beta_003125to16.α) + 1)
data8[1, 1] = "w(ω₀)|T(ħω₀/kB)↓|α→"
data8[2:end, 1] = ustrip(1 ./ frohlich_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data8[1, 2:end] = 0.1:0.1:12
data8[2:end, 2:end] = reverse(ustrip(frohlich_3d_alpha_0to12_beta_003125to16.w[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-w-alpha-0to12-temp-00625to32.dat", data8)
data8
end

# ╔═╡ 59947a9f-4599-4af8-a9c5-1f5f249ddc45
md"""
## M = v²/w² m₀
"""

# ╔═╡ 1b6463e2-d744-4d7a-95fe-9d09ecd65e7f
md"""
### 2D
"""

# ╔═╡ 44c47eaa-bef8-4c50-bacd-b10866558caa
begin
data9 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data9[1, 1] = "M(m₀)|T(ħω₀/kB)↓|α→"
data9[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data9[1, 2:end] = 0.1:0.1:12
data9[2:end, 2:end] = reverse(ustrip(polaron_mass(frohlich_2d_alpha_0to12_beta_003125to16)[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-mass-alpha-0to12-temp-00625to32.dat", data9)
data9
end

# ╔═╡ 90a2c010-0650-4952-96d7-f8b03b7e8221
md"""
### 3D
"""

# ╔═╡ 9f1127a5-e401-4630-a429-71a862a54b44
begin
data10 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_3d_alpha_0to12_beta_003125to16.α) + 1)
data10[1, 1] = "M(m₀)|T(ħω₀/kB)↓|α→"
data10[2:end, 1] = ustrip(1 ./ frohlich_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data10[1, 2:end] = 0.1:0.1:12
data10[2:end, 2:end] = reverse(ustrip(polaron_mass(frohlich_3d_alpha_0to12_beta_003125to16)[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-mass-alpha-0to12-temp-00625to32.dat", data10)
data10
end

# ╔═╡ 367c3815-0a85-493e-93bf-1d567ed46fe0
md"""
## κ = v²-w² m₀ω₀²
"""

# ╔═╡ edb19919-7952-4578-93dd-78c1373d8036
md"""
### 2D
"""

# ╔═╡ a54fbca9-2752-45d1-9ded-06218ad89afe
begin
data11 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data11[1, 1] = "κ(m₀ω₀²)|T(ħω₀/kB)↓|α→"
data11[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data11[1, 2:end] = 0.1:0.1:12
data11[2:end, 2:end] = reverse(ustrip(polaron_spring(frohlich_2d_alpha_0to12_beta_003125to16)[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-spring-alpha-0to12-temp-00625to32.dat", data11)
data11
end

# ╔═╡ 52ce12a6-e895-4059-87ab-a1a26c6e05ea
md"""
### 3D
"""

# ╔═╡ fa667323-3183-41d0-8958-1ff5cefe7e9f
begin
data12 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data12[1, 1] = "κ(m₀ω₀²)|T(ħω₀/kB)↓|α→"
data12[2:end, 1] = ustrip(1 ./ frohlich_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data12[1, 2:end] = 0.1:0.1:12
data12[2:end, 2:end] = reverse(ustrip(polaron_spring(frohlich_3d_alpha_0to12_beta_003125to16)[:,1,:])',dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-spring-alpha-0to12-temp-00625to32.dat", data12)
data12
end

# ╔═╡ 76d0276b-47fd-4fb9-ae06-21c3d64a5f8a
md"""
## R = D/2 v/(v²-w²) coth(ħβv/2) r₀=√(ħ/m₀ω₀)
"""

# ╔═╡ 28e72ed4-aa97-482c-8dd4-35add756f650
md"""
### 2D
"""

# ╔═╡ e2b40339-f0da-4ee3-8017-1e18708b9175
begin
data13 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data13[1, 1] = "R(r₀)|T(ħω₀/kB)↓|α→"
data13[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data13[1, 2:end] = 0.1:0.1:12
data13[2:end, 2:end] = reverse(ustrip(polaron_radius2(frohlich_2d_alpha_0to12_beta_003125to16,dims=2)),dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-radius-alpha-0to12-temp-00625to32.dat", data13)
data13
end

# ╔═╡ f5fb95ae-1be6-42cc-9c9f-8343590a1a4c
md"""
### 3D
"""

# ╔═╡ 658e98cb-5409-4b4b-9299-5d1ba8247cc1
begin
data14 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_3d_alpha_0to12_beta_003125to16.α) + 1)
data14[1, 1] = "R(r₀)|T(ħω₀/kB)↓|α→"
data14[2:end, 1] = ustrip(1 ./ frohlich_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data14[1, 2:end] = 0.1:0.1:12
data14[2:end, 2:end] = reverse(ustrip(polaron_radius2(frohlich_3d_alpha_0to12_beta_003125to16)),dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-radius-alpha-0to12-temp-00625to32.dat", data14)
data14
end

# ╔═╡ 790a9562-ce93-497b-9705-8fdb8985224c
md"""
## μ = 1/ℑΣ eω₀/m₀
"""

# ╔═╡ d6aa74e2-d5d7-46d3-a7ca-27bba134f2d0
md"""
### 2D
"""

# ╔═╡ b3cca528-a7cb-442d-924e-d2fe83fad373
begin
data15 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_2d_alpha_0to12_beta_003125to16.α) + 1)
data15[1, 1] = "μ(eω₀/m₀)|T(ħω₀/kB)↓|α→"
data15[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data15[1, 2:end] = 0.1:0.1:12
data15[2:end, 2:end] = reverse(ustrip(polaron_mobility(frohlich_2d_alpha_0to12_beta_003125to16))',dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat", data15)
data15
end

# ╔═╡ 3f1ad58b-71a0-4d65-8139-b79d37d765e0
md"""
### 3D
"""

# ╔═╡ 1f4a385a-3661-4bba-9d57-1615a0e4d60b
begin
data16 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_003125to16.β) + 1, length(frohlich_3d_alpha_0to12_beta_003125to16.α) + 1)
data16[1, 1] = "μ(eω₀/m₀)|T(ħω₀/kB)↓|α→"
data16[2:end, 1] = ustrip(1 ./ frohlich_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data16[1, 2:end] = 0.1:0.1:12
data16[2:end, 2:end] = reverse(ustrip(polaron_mobility(frohlich_3d_alpha_0to12_beta_003125to16))',dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-mobility-alpha-0to12-temp-00625to32.dat", data16)
data16
end

# ╔═╡ 99354529-397a-41da-9bd0-c328a02025e1
md"""
# 2D and 3D Frohlich α = 0.05:0.05:12, T = 0.48 ħω₀/kB, Ω = 0.01:0.01:30.0 ω₀
"""

# ╔═╡ 0636ff35-bfea-4b0e-9e3b-af3559b1468c
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 30734s
begin
frohlich_2d_alpha_0to12_beta_100_freq_0to30 = frohlich(LinRange(12,0.05,240), 1.0, 100.0, 0.01:0.01:30; dims = 2)
save_frohlich(frohlich_2d_alpha_0to12_beta_100_freq_0to30, "data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-100-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 4bfa6a60-757d-459b-8cbc-3de399639d35
frohlich_2d_alpha_0to12_beta_100_freq_0to30 = load_frohlich("data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-100-freq-0to30.jld")

# ╔═╡ 99859951-fce1-4224-9961-cb8bfb3a0ddb
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 8032s
begin
frohlich_3d_alpha_0to12_beta_100_freq_0to30 = frohlich(LinRange(12,0.05,240), 1.0, 100.0, 0.01:0.01:30)
save_frohlich(frohlich_3d_alpha_0to12_beta_100_freq_0to30, "data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-100-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 9ae9e209-3db0-46e9-8e36-b4394a4d9417
frohlich_3d_alpha_0to12_beta_100_freq_0to30 = load_frohlich("data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-100-freq-0to30.jld")

# ╔═╡ 724fa1d6-cf9a-4c41-bff1-fc9108086d14
md"""
## Σ ω₀
"""

# ╔═╡ 5bf8f537-acac-4a5b-9fc3-f453d2ec2e29
md"""
### ℜΣ
"""

# ╔═╡ c0f61a12-276a-451f-a94c-c52855af8d04
begin
data190 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data190[1, 1] = "ℜΣ(ω₀)|Ω(ω₀)↓|α→"
data190[2:end, 1] = 0.01:0.01:30.0
data190[1, 2:end] = 0.05:0.05:12
data190[2:end, 2:end] = reverse(ustrip(real.(frohlich_2d_alpha_0to12_beta_100_freq_0to30.Σ)'),dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-real-memory-alpha-0to12-beta-100-freq-0to30.dat", data190)
data190
end

# ╔═╡ 6449a1f8-4faf-4216-80a1-9f2b63ca5d2e
begin
data19 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data19[1, 1] = "ℜΣ(ω₀)|Ω(ω₀)↓|α→"
data19[2:end, 1] = 0.01:0.01:30.0
data19[1, 2:end] = 0.05:0.05:12
data19[2:end, 2:end] = reverse(ustrip(real.(frohlich_3d_alpha_0to12_beta_100_freq_0to30.Σ)'),dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-real-memory-alpha-0to12-beta-100-freq-0to30.dat", data19)
data19
end

# ╔═╡ d843c9af-7833-40ae-9391-90b3c55047ad
md"""
### ℑΣ
"""

# ╔═╡ 61ae1a96-cb58-4e0c-af35-98d257fa5bbe
begin
data200 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data200[1, 1] = "ℑΣ(ω₀)|Ω(ω₀)↓|α→"
data200[2:end, 1] = 0.01:0.01:30.0
data200[1, 2:end] = 0.05:0.05:12
data200[2:end, 2:end] = reverse(ustrip(abs.(imag.(frohlich_2d_alpha_0to12_beta_100_freq_0to30.Σ)')),dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat", data200)
data200
end

# ╔═╡ f20bbb33-995c-46a7-871f-9090246a14cc
begin
data20 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data20[1, 1] = "ℑΣ(ω₀)|Ω(ω₀)↓|α→"
data20[2:end, 1] = 0.01:0.01:30.0
data20[1, 2:end] = 0.05:0.05:12
data20[2:end, 2:end] = reverse(ustrip(abs.(imag.(frohlich_3d_alpha_0to12_beta_100_freq_0to30.Σ)')),dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat", data20)
data20
end

# ╔═╡ c7573414-d8b7-47db-8f01-9b3d7728117f
md"""
## σ = -im/(Ω+Σ) e²ω₀/m₀r₀²=e²ω₀²/ħ
"""

# ╔═╡ 9c8ecb33-095b-4526-b38f-f7ef7f3a50e9
md"""
### ℜσ
"""

# ╔═╡ 17200fcc-539c-48ec-880c-a276fc593415
begin
data170 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data170[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|Ω(ω₀)↓|α→"
data170[2:end, 1] = 0.01:0.01:30.0
data170[1, 2:end] = 0.05:0.05:12
data170[2:end, 2:end] = reverse(ustrip(abs.(real.(polaron_conductivity(frohlich_2d_alpha_0to12_beta_100_freq_0to30)))),dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-real-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data170)
data170
end

# ╔═╡ 2f816ee0-d6a2-4eb3-9027-a54461dd83bd
begin
data17 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data17[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|Ω(ω₀)↓|α→"
data17[2:end, 1] = 0.01:0.01:30.0
data17[1, 2:end] = 0.05:0.05:12
data17[2:end, 2:end] = reverse(ustrip(abs.(real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_freq_0to30)))),dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-real-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data17)
data17
end

# ╔═╡ 962d8925-2b86-4af1-8971-21e9832c5faf
md"""
### ℑσ
"""

# ╔═╡ d2682bc3-1d04-4619-b69d-c713852b7e11
begin
data180 = Matrix{Any}(undef, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_2d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data180[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|Ω(ω₀)↓|α→"
data180[2:end, 1] = 0.01:0.01:30.0
data180[1, 2:end] = 0.05:0.05:12
data180[2:end, 2:end] = reverse(ustrip(imag.(polaron_conductivity(frohlich_2d_alpha_0to12_beta_100_freq_0to30))),dims=2)
writedlm("data/frohlich/variational/model/frohlich-2d-imag-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data180)
data180
end

# ╔═╡ e8a3b9af-24b0-4b6c-974f-5049ae1e4a43
begin
data18 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_100_freq_0to30.α) + 1)
data18[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|Ω(ω₀)↓|α→"
data18[2:end, 1] = 0.01:0.01:30.0
data18[1, 2:end] = 0.05:0.05:12
data18[2:end, 2:end] = reverse(ustrip(imag.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_freq_0to30))),dims=2)
writedlm("data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data18)
data18
end

# ╔═╡ d83ae509-e1cd-4178-94bd-4037a68f2ac8
md"""
# 3D Frohlich α = 1.0:7.0, T = LogRange(0.0625,32.0,128) ħω₀/kB, Ω = 0.01:0.01:30.0 ω₀
"""

# ╔═╡ e1311bc1-7cb5-41be-8a1f-5104a8630500
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 1816s
for α in 1:5
	@eval $(Symbol("frohlich_2d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(frohlich(α, 1.0, LogRange(0.0325,16.0,128), 0.01:0.01:30.0; dims = 2))
	save_frohlich(eval(Symbol("frohlich_2d_alpha_$(α)_beta_00325to16_freq_0to30")), "data/frohlich/variational/model/frohlich-2d-alpha-$(α)-beta-00325to16-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ bb1efbcf-7728-4f18-9a5c-d73315c40459
for α in 1:5
	@eval $(Symbol("frohlich_2d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(load_frohlich("data/frohlich/variational/model/frohlich-2d-alpha-$α-beta-00325to16-freq-0to30.jld"))
end

# ╔═╡ 4e4b75c3-582b-45b0-b389-49d9e2552817
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 3725s
for α in 1:7
	@eval $(Symbol("frohlich_3d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(frohlich(α, 1.0, LogRange(0.0325,16.0,128), 0.01:0.01:30.0))
	save_frohlich(eval(Symbol("frohlich_3d_alpha_$(α)_beta_00325to16_freq_0to30")), "data/frohlich/variational/model/frohlich-3d-alpha-$(α)-beta-00325to16-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 1d66f7a6-f8fa-478f-b007-956f24d047f7
for α in 1:7
	@eval $(Symbol("frohlich_3d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(load_frohlich("data/frohlich/variational/model/frohlich-3d-alpha-$α-beta-00325to16-freq-0to30.jld"))
end

# ╔═╡ 9804abff-66b2-40af-8c16-a099e7600b73
md"""
## Σ ω₀
"""

# ╔═╡ 5fef1d73-a691-438b-8634-293c098f1522
md"""
### ℜΣ
"""

# ╔═╡ 092ffc9d-e65a-4ba3-ad7b-ee346f6ff06d
for α in 1:5
	f = eval(Symbol("frohlich_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data22 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data22[1, 1] = "ℜΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data22[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data22[1, 2:end] = 0.01:0.01:30
	data22[2:end, 2:end] = ustrip(real.(f.Σ))
	writedlm("data/frohlich/variational/model/frohlich-2d-real-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data22)
end

# ╔═╡ 9ce6fad2-9746-4605-b967-5447f8d0654d
for α in 1:7
	f = eval(Symbol("frohlich_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data22 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data22[1, 1] = "ℜΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data22[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data22[1, 2:end] = 0.01:0.01:30
	data22[2:end, 2:end] = ustrip(real.(f.Σ))
	writedlm("data/frohlich/variational/model/frohlich-3d-real-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data22)
end

# ╔═╡ 45377ee4-b055-4ed4-8bb3-2b7321c4e709
md"""
### ℑΣ
"""

# ╔═╡ 5b5a93ac-13c7-49de-bf6d-e5e546c9c2fb
for α in 1:5
	f = eval(Symbol("frohlich_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data21 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data21[1, 1] = "ℑΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data21[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data21[1, 2:end] = 0.01:0.01:30
	data21[2:end, 2:end] = ustrip(abs.(imag.(f.Σ)))
	writedlm("data/frohlich/variational/model/frohlich-2d-imag-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data21)
end

# ╔═╡ 20f9700f-3cb9-4197-8b9a-b9fe71cc748c
for α in 1:7
	f = eval(Symbol("frohlich_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data21 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data21[1, 1] = "ℑΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data21[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data21[1, 2:end] = 0.01:0.01:30
	data21[2:end, 2:end] = ustrip(abs.(imag.(f.Σ)))
	writedlm("data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data21)
end

# ╔═╡ b17f1e29-cc52-455b-a6ec-aec4dc70a372
md"""
## σ = -im/(Ω+Σ) e²ω₀/m₀r₀²=e²ω₀²/ħ
"""

# ╔═╡ 4d3645bb-cdf3-4e2c-9edd-c6c5be475107
md"""
### ℜσ
"""

# ╔═╡ 7168a55f-e46d-4e9a-bd5c-e0daf312d126
for α in 1:5
	f = eval(Symbol("frohlich_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data23 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data23[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data23[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data23[1, 2:end] = 0.01:0.01:30
	data23[2:end, 2:end] = ustrip(abs.(real.(polaron_conductivity(f))))'
	writedlm("data/frohlich/variational/model/frohlich-2d-real-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data23)
end

# ╔═╡ 84f94766-918f-40d6-a51a-2124f9e0df86
for α in 1:7
	f = eval(Symbol("frohlich_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data23 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data23[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data23[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data23[1, 2:end] = 0.01:0.01:30
	data23[2:end, 2:end] = ustrip(abs.(real.(polaron_conductivity(f))))'
	writedlm("data/frohlich/variational/model/frohlich-3d-real-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data23)
end

# ╔═╡ fdc5be06-e1d0-4346-9402-0428fa451f8f
md"""
### ℑσ
"""

# ╔═╡ 829640a1-7787-4c46-8954-13c96c7bfe29
for α in 1:5
	f = eval(Symbol("frohlich_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data24 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data24[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data24[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data24[1, 2:end] = 0.01:0.01:30
	data24[2:end, 2:end] = ustrip(imag.(polaron_conductivity(f)))'
	writedlm("data/frohlich/variational/model/frohlich-2d-imag-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data24)
end

# ╔═╡ 7f1fcea4-1933-442b-bdab-88c37e5b79b0
for α in 1:7
	f = eval(Symbol("frohlich_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data24 = Matrix{Any}(undef, length(f.β) + 1, length(f.Ω) + 1)
	data24[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data24[2:end, 1] = ustrip(1 ./ f.β ./ kB_pu)
	data24[1, 2:end] = 0.01:0.01:30
	data24[2:end, 2:end] = ustrip(imag.(polaron_conductivity(f)))'
	writedlm("data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data24)
end

# ╔═╡ Cell order:
# ╟─bfd6e85c-48e2-11ef-2751-3147edce85d8
# ╟─68c1b156-ee71-4717-94d6-914500fab1ff
# ╟─52e5925c-a82a-4799-97b4-a96a0798fa9b
# ╟─b12c5fa7-a588-4ce1-a4b5-1cf598380903
# ╟─4d74ec93-b41d-40af-b8e6-0652b2f10faf
# ╟─83427388-3d3f-4f05-867a-6984454cfdf9
# ╟─8564e569-76ed-439a-aa13-291538cfff13
# ╟─828f01c9-5726-4c8a-a5b1-bd6e83f29662
# ╟─6088c544-f44b-49c1-a503-ca881418e2f9
# ╟─8420455b-5772-413d-8e19-4e541e379142
# ╟─e03e239c-a77d-4591-aaa0-e02ce4cc0292
# ╟─9d95d5b9-1eee-4d7c-acd2-a2177f620fb7
# ╟─174468b5-1853-4eeb-858f-3a65f9c8a822
# ╟─c19ed360-aa81-4ed1-a5b6-0359b01e9a5f
# ╟─2cb035de-37aa-434b-b61e-a9106bac6b86
# ╟─90a6cf36-3b99-4841-bb6a-5c0192c11c37
# ╟─05caea00-de80-4357-83b3-8d884b8e0573
# ╟─9a91d18a-b49f-46b5-96e1-f577b142ceac
# ╟─477eaae2-731a-439a-9888-1d00f42943a8
# ╟─92b52583-b5db-43cd-8ee9-06098deece54
# ╟─4589e7b4-5bfe-437f-927b-35626eb52b84
# ╟─ae307e8c-13c4-4d5b-8f99-4c7242df1565
# ╟─808e75f2-7956-4eef-bec8-0564c8750aec
# ╟─a9de638c-c8f7-438a-89fc-e362d7e725cb
# ╟─8810d4e6-4ce8-4c29-af10-e34d7470a062
# ╟─6a881a63-e422-4389-b25b-23a2c6d82d29
# ╟─c8d1b9fc-9c07-4ebb-8694-05f9586943eb
# ╟─f134f1d4-33c1-48c8-89ec-0f9f2d457e3e
# ╟─99faadff-2dd4-46ce-9a8c-d8faf850c9b7
# ╟─958e1299-defb-410d-ac85-c176f621462c
# ╟─caa54305-0f13-47a4-850d-a3b14d975ad1
# ╟─e5d9f83f-1c5d-4e0a-9842-554846c227dc
# ╟─13b93723-dcd3-4bc5-9935-30635ec2bbab
# ╟─1e306f70-be2b-4704-aa70-9d94d85d8859
# ╟─c417194c-084c-440e-9bed-a5d7b9219ac8
# ╟─2f330d5d-ebd5-4e56-b9d6-5e4dc8c48ac4
# ╟─670e50ac-afea-4e62-aef7-9bd12a8541d5
# ╟─e8c3adce-01b9-4ffe-a946-585cce2e4ae4
# ╟─c3a3a96a-cc7c-4f13-9a01-14473007fcd9
# ╟─59947a9f-4599-4af8-a9c5-1f5f249ddc45
# ╟─1b6463e2-d744-4d7a-95fe-9d09ecd65e7f
# ╟─44c47eaa-bef8-4c50-bacd-b10866558caa
# ╟─90a2c010-0650-4952-96d7-f8b03b7e8221
# ╟─9f1127a5-e401-4630-a429-71a862a54b44
# ╟─367c3815-0a85-493e-93bf-1d567ed46fe0
# ╟─edb19919-7952-4578-93dd-78c1373d8036
# ╟─a54fbca9-2752-45d1-9ded-06218ad89afe
# ╟─52ce12a6-e895-4059-87ab-a1a26c6e05ea
# ╟─fa667323-3183-41d0-8958-1ff5cefe7e9f
# ╟─76d0276b-47fd-4fb9-ae06-21c3d64a5f8a
# ╟─28e72ed4-aa97-482c-8dd4-35add756f650
# ╟─e2b40339-f0da-4ee3-8017-1e18708b9175
# ╟─f5fb95ae-1be6-42cc-9c9f-8343590a1a4c
# ╟─658e98cb-5409-4b4b-9299-5d1ba8247cc1
# ╟─790a9562-ce93-497b-9705-8fdb8985224c
# ╟─d6aa74e2-d5d7-46d3-a7ca-27bba134f2d0
# ╟─b3cca528-a7cb-442d-924e-d2fe83fad373
# ╟─3f1ad58b-71a0-4d65-8139-b79d37d765e0
# ╟─1f4a385a-3661-4bba-9d57-1615a0e4d60b
# ╟─99354529-397a-41da-9bd0-c328a02025e1
# ╟─0636ff35-bfea-4b0e-9e3b-af3559b1468c
# ╟─4bfa6a60-757d-459b-8cbc-3de399639d35
# ╟─99859951-fce1-4224-9961-cb8bfb3a0ddb
# ╟─9ae9e209-3db0-46e9-8e36-b4394a4d9417
# ╟─724fa1d6-cf9a-4c41-bff1-fc9108086d14
# ╟─5bf8f537-acac-4a5b-9fc3-f453d2ec2e29
# ╟─c0f61a12-276a-451f-a94c-c52855af8d04
# ╟─6449a1f8-4faf-4216-80a1-9f2b63ca5d2e
# ╟─d843c9af-7833-40ae-9391-90b3c55047ad
# ╟─61ae1a96-cb58-4e0c-af35-98d257fa5bbe
# ╟─f20bbb33-995c-46a7-871f-9090246a14cc
# ╟─c7573414-d8b7-47db-8f01-9b3d7728117f
# ╟─9c8ecb33-095b-4526-b38f-f7ef7f3a50e9
# ╟─17200fcc-539c-48ec-880c-a276fc593415
# ╟─2f816ee0-d6a2-4eb3-9027-a54461dd83bd
# ╟─962d8925-2b86-4af1-8971-21e9832c5faf
# ╟─d2682bc3-1d04-4619-b69d-c713852b7e11
# ╟─e8a3b9af-24b0-4b6c-974f-5049ae1e4a43
# ╟─d83ae509-e1cd-4178-94bd-4037a68f2ac8
# ╠═e1311bc1-7cb5-41be-8a1f-5104a8630500
# ╠═bb1efbcf-7728-4f18-9a5c-d73315c40459
# ╠═4e4b75c3-582b-45b0-b389-49d9e2552817
# ╟─1d66f7a6-f8fa-478f-b007-956f24d047f7
# ╟─9804abff-66b2-40af-8c16-a099e7600b73
# ╟─5fef1d73-a691-438b-8634-293c098f1522
# ╠═092ffc9d-e65a-4ba3-ad7b-ee346f6ff06d
# ╠═9ce6fad2-9746-4605-b967-5447f8d0654d
# ╟─45377ee4-b055-4ed4-8bb3-2b7321c4e709
# ╠═5b5a93ac-13c7-49de-bf6d-e5e546c9c2fb
# ╠═20f9700f-3cb9-4197-8b9a-b9fe71cc748c
# ╟─b17f1e29-cc52-455b-a6ec-aec4dc70a372
# ╟─4d3645bb-cdf3-4e2c-9edd-c6c5be475107
# ╠═7168a55f-e46d-4e9a-bd5c-e0daf312d126
# ╠═84f94766-918f-40d6-a51a-2124f9e0df86
# ╟─fdc5be06-e1d0-4346-9402-0428fa451f8f
# ╠═829640a1-7787-4c46-8954-13c96c7bfe29
# ╠═7f1fcea4-1933-442b-bdab-88c37e5b79b0
