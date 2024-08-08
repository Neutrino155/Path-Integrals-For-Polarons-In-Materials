### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ bfd6e85c-48e2-11ef-2751-3147edce85d8
using Pkg, Unitful, Plots, DelimitedFiles

# ╔═╡ 35632de9-363d-459e-87b0-4bfaaccddb1f
Pkg.activate(".")

# ╔═╡ 52e5925c-a82a-4799-97b4-a96a0798fa9b
using PolaronThesis

# ╔═╡ b12c5fa7-a588-4ce1-a4b5-1cf598380903
md"""
# Polaron Equations
"""

# ╔═╡ 4d74ec93-b41d-40af-b8e6-0652b2f10faf
polaron_mass(holstein) = holstein.v.^2 ./ holstein.w.^2

# ╔═╡ 83427388-3d3f-4f05-867a-6984454cfdf9
polaron_spring(holstein) = holstein.v.^2 .- holstein.w.^2

# ╔═╡ 8564e569-76ed-439a-aa13-291538cfff13
polaron_radius(holstein; dims=3) = dims / 2 .* holstein.v ./ (holstein.v.^2 .- holstein.w.^2) .* coth.(pustrip.(holstein.β .* holstein.v ./ 2))

# ╔═╡ 828f01c9-5726-4c8a-a5b1-bd6e83f29662
polaron_radius2(holstein; dims=3) = dims / 2 .* pustrip.(holstein.v[:,1,1,:])' ./ (pustrip.(holstein.v[:,1,1,:])'.^2 .- pustrip.(holstein.w[:,1,1,:])'.^2) .* coth.(pustrip.(holstein.β .* pustrip.(holstein.v[:,1,1,:])' ./ 2))

# ╔═╡ 6088c544-f44b-49c1-a503-ca881418e2f9
polaron_mobility(holstein) = abs.(1 ./ imag.(holstein.Σ))

# ╔═╡ 8420455b-5772-413d-8e19-4e541e379142
polaron_conductivity(holstein) = -im ./ (holstein.Ω .+ holstein.Σ')

# ╔═╡ e03e239c-a77d-4591-aaa0-e02ce4cc0292
md"""
# 1D, 2D and 3D Holstein α = 0.1:0.1:12, ω = 0.1:0.1:2.0
"""

# ╔═╡ 6c2184de-0a93-43df-8b78-778b23844094
md"""
## 1D
"""

# ╔═╡ 4ee11b9d-f78c-437d-8e07-0c1fae31c2d5
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 14s
begin
holstein_1d_alpha_0to12_beta_inf = holstein(12.0:-0.1:0.1, 0.1:0.1:2.0, 1.0, dims = 1)
save_holstein(holstein_1d_alpha_0to12_beta_inf, "data/holstein/variational/model/holstein-1d-alpha-0to12-beta-inf")
end
  ╠═╡ =#

# ╔═╡ 59777dfa-61f6-49fe-8b7b-5a68dd691c86
holstein_1d_alpha_0to12_beta_inf = load_holstein("data/holstein/variational/model/holstein-1d-alpha-0to12-beta-inf.jld")

# ╔═╡ 057fde9d-dffc-4ca0-acde-813ba4c8966a
begin
data00 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_inf.α) + 1,  length(holstein_1d_alpha_0to12_beta_inf.ω) + 1)
data00[1, 1] = "E(ħω₀)|α↓|ω→"
data00[1, 2:end] = holstein_1d_alpha_0to12_beta_inf.ω
data00[2:end, 1] = 0.1:0.1:12
data00[2:end, 2:end] = reverse(ustrip(holstein_1d_alpha_0to12_beta_inf.E[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-energy-alpha-0to12-beta-inf.dat", data00)
data00
end

# ╔═╡ 45f46b37-e45d-4a01-ba51-9f00329876e6
begin
data01 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_inf.α) + 1,  length(holstein_1d_alpha_0to12_beta_inf.ω) + 1)
data01[1, 1] = "v(ω₀)|α↓|ω→"
data01[1, 2:end] = holstein_1d_alpha_0to12_beta_inf.ω
data01[2:end, 1] = 0.1:0.1:12
data01[2:end, 2:end] = reverse(ustrip(holstein_1d_alpha_0to12_beta_inf.v[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-v-alpha-0to12-beta-inf.dat", data01)
data01
end

# ╔═╡ 7b032f37-8a76-4db8-82c4-e7bf6719a031
begin
data02 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_inf.α) + 1,  length(holstein_1d_alpha_0to12_beta_inf.ω) + 1)
data02[1, 1] = "w(ω₀)|α↓|ω→"
data02[1, 2:end] = holstein_1d_alpha_0to12_beta_inf.ω
data02[2:end, 1] = 0.1:0.1:12
data02[2:end, 2:end] = reverse(ustrip(holstein_1d_alpha_0to12_beta_inf.w[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-w-alpha-0to12-beta-inf.dat", data02)
data02
end

# ╔═╡ 63b234ce-6b5f-4645-a3bb-04225f9575e7
begin
data03 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_inf.α) + 1,  length(holstein_1d_alpha_0to12_beta_inf.ω) + 1)
data03[1, 1] = "M(m₀)|α↓|ω→"
data03[1, 2:end] = holstein_1d_alpha_0to12_beta_inf.ω
data03[2:end, 1] = 0.1:0.1:12
data03[2:end, 2:end] = reverse(ustrip(polaron_mass(holstein_1d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-mass-alpha-0to12-beta-inf.dat", data03)
data03
end

# ╔═╡ 8cd17acc-9a62-45d6-8016-1af47643403a
begin
data04 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_inf.α) + 1,  length(holstein_1d_alpha_0to12_beta_inf.ω) + 1)
data04[1, 1] = "κ(m₀ω₀²)|α↓|ω→"
data04[1, 2:end] = holstein_1d_alpha_0to12_beta_inf.ω
data04[2:end, 1] = 0.1:0.1:12
data04[2:end, 2:end] = reverse(ustrip(polaron_spring(holstein_1d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-spring-alpha-0to12-beta-inf.dat", data04)
data04
end

# ╔═╡ dd2c5cfd-b530-4380-9b5b-29739dc7217e
begin
data05 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_inf.α) + 1,  length(holstein_1d_alpha_0to12_beta_inf.ω) + 1)
data05[1, 1] = "R(r₀)|α↓|ω→"
data05[1, 2:end] = holstein_1d_alpha_0to12_beta_inf.ω
data05[2:end, 1] = 0.1:0.1:12
data05[2:end, 2:end] = reverse(ustrip(polaron_radius(holstein_1d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-radius-alpha-0to12-beta-inf.dat", data05)
data05
end

# ╔═╡ 9d95d5b9-1eee-4d7c-acd2-a2177f620fb7
md"""
## 2D
"""

# ╔═╡ 174468b5-1853-4eeb-858f-3a65f9c8a822
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 22s
begin
holstein_2d_alpha_0to12_beta_inf = holstein(12.0:-0.1:0.1, 0.1:0.1:2.0, 1.0; dims = 2)
save_holstein(holstein_2d_alpha_0to12_beta_inf, "data/holstein/variational/model/holstein-2d-alpha-0to12-beta-inf")
end
  ╠═╡ =#

# ╔═╡ c19ed360-aa81-4ed1-a5b6-0359b01e9a5f
holstein_2d_alpha_0to12_beta_inf = load_holstein("data/holstein/variational/model/holstein-2d-alpha-0to12-beta-inf.jld")

# ╔═╡ 2cb035de-37aa-434b-b61e-a9106bac6b86
begin
data06 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_inf.α) + 1,  length(holstein_2d_alpha_0to12_beta_inf.ω) + 1)
data06[1, 1] = "E(ħω₀)|α↓|ω→"
data06[1, 2:end] = holstein_2d_alpha_0to12_beta_inf.ω
data06[2:end, 1] = 0.1:0.1:12
data06[2:end, 2:end] = reverse(ustrip(holstein_2d_alpha_0to12_beta_inf.E[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-energy-alpha-0to12-beta-inf.dat", data06)
data06
end

# ╔═╡ 7aec69dd-b635-4ae8-9be2-92f0fa8a6851
begin
data07 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_inf.α) + 1,  length(holstein_2d_alpha_0to12_beta_inf.ω) + 1)
data07[1, 1] = "v(ω₀)|α↓|ω→"
data07[1, 2:end] = holstein_2d_alpha_0to12_beta_inf.ω
data07[2:end, 1] = 0.1:0.1:12
data07[2:end, 2:end] = reverse(ustrip(holstein_2d_alpha_0to12_beta_inf.v[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-v-alpha-0to12-beta-inf.dat", data07)
data07
end

# ╔═╡ be95e280-e5bd-4669-85f5-dbcb835b7b51
begin
data08 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_inf.α) + 1,  length(holstein_2d_alpha_0to12_beta_inf.ω) + 1)
data08[1, 1] = "w(ω₀)|α↓|ω→"
data08[1, 2:end] = holstein_2d_alpha_0to12_beta_inf.ω
data08[2:end, 1] = 0.1:0.1:12
data08[2:end, 2:end] = reverse(ustrip(holstein_2d_alpha_0to12_beta_inf.w[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-w-alpha-0to12-beta-inf.dat", data08)
data08
end

# ╔═╡ 82f79186-360b-4b78-a709-407bde56621c
begin
data09 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_inf.α) + 1,  length(holstein_2d_alpha_0to12_beta_inf.ω) + 1)
data09[1, 1] = "M(m₀)|α↓|ω→"
data09[1, 2:end] = holstein_2d_alpha_0to12_beta_inf.ω
data09[2:end, 1] = 0.1:0.1:12
data09[2:end, 2:end] = reverse(ustrip(polaron_mass(holstein_2d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-mass-alpha-0to12-beta-inf.dat", data09)
data09
end

# ╔═╡ 337a233b-3707-4103-b555-2e4fa94a1d14
begin
data010 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_inf.α) + 1,  length(holstein_2d_alpha_0to12_beta_inf.ω) + 1)
data010[1, 1] = "κ(m₀ω₀²)|α↓|ω→"
data010[1, 2:end] = holstein_2d_alpha_0to12_beta_inf.ω
data010[2:end, 1] = 0.1:0.1:12
data010[2:end, 2:end] = reverse(ustrip(polaron_spring(holstein_2d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-spring-alpha-0to12-beta-inf.dat", data010)
data010
end

# ╔═╡ 1f5dc90e-e61c-4213-8c02-09030b011a44
begin
data011 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_inf.α) + 1,  length(holstein_2d_alpha_0to12_beta_inf.ω) + 1)
data011[1, 1] = "R(r₀)|α↓|ω→"
data011[1, 2:end] = holstein_2d_alpha_0to12_beta_inf.ω
data011[2:end, 1] = 0.1:0.1:12
data011[2:end, 2:end] = reverse(ustrip(polaron_radius(holstein_2d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-adius-alpha-0to12-beta-inf.dat", data011)
data011
end

# ╔═╡ 90a6cf36-3b99-4841-bb6a-5c0192c11c37
md"""
## 3D
"""

# ╔═╡ 05caea00-de80-4357-83b3-8d884b8e0573
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 23s
begin
holstein_3d_alpha_0to12_beta_inf = holstein(12.0:-0.1:0.1, 0.1:0.1:2.0, 1.0; dims = 3)
save_holstein(holstein_3d_alpha_0to12_beta_inf, "data/holstein/variational/model/holstein-3d-alpha-0to12-beta-inf")
end
  ╠═╡ =#

# ╔═╡ 9a91d18a-b49f-46b5-96e1-f577b142ceac
holstein_3d_alpha_0to12_beta_inf = load_holstein("data/holstein/variational/model/holstein-3d-alpha-0to12-beta-inf.jld")

# ╔═╡ 477eaae2-731a-439a-9888-1d00f42943a8
begin
data012 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_inf.α) + 1,  length(holstein_3d_alpha_0to12_beta_inf.ω) + 1)
data012[1, 1] = "E(ħω₀)|α↓|ω→"
data012[1, 2:end] = holstein_3d_alpha_0to12_beta_inf.ω
data012[2:end, 1] = 0.1:0.1:12
data012[2:end, 2:end] = reverse(ustrip(holstein_3d_alpha_0to12_beta_inf.E[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-energy-alpha-0to12-beta-inf.dat", data012)
data012
end

# ╔═╡ 03d5c9cb-fe71-44c9-9838-f95283a3f874
begin
data013 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_inf.α) + 1,  length(holstein_3d_alpha_0to12_beta_inf.ω) + 1)
data013[1, 1] = "v(ω₀)|α↓|ω→"
data013[1, 2:end] = holstein_3d_alpha_0to12_beta_inf.ω
data013[2:end, 1] = 0.1:0.1:12
data013[2:end, 2:end] = reverse(ustrip(holstein_3d_alpha_0to12_beta_inf.v[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-v-alpha-0to12-beta-inf.dat", data013)
data013
end

# ╔═╡ d04fe83b-a165-4c43-8316-906c1944fe5c
begin
data014 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_inf.α) + 1,  length(holstein_3d_alpha_0to12_beta_inf.ω) + 1)
data014[1, 1] = "w(ω₀)|α↓|ω→"
data014[1, 2:end] = holstein_3d_alpha_0to12_beta_inf.ω
data014[2:end, 1] = 0.1:0.1:12
data014[2:end, 2:end] = reverse(ustrip(holstein_3d_alpha_0to12_beta_inf.w[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-w-alpha-0to12-beta-inf.dat", data014)
data014
end

# ╔═╡ 0d60752c-f34a-4031-a1f4-57b8aa8c1717
begin
data015 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_inf.α) + 1,  length(holstein_3d_alpha_0to12_beta_inf.ω) + 1)
data015[1, 1] = "M(m₀)|α↓|ω→"
data015[1, 2:end] = holstein_3d_alpha_0to12_beta_inf.ω
data015[2:end, 1] = 0.1:0.1:12
data015[2:end, 2:end] = reverse(ustrip(polaron_mass(holstein_3d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-mass-alpha-0to12-beta-inf.dat", data015)
data015
end

# ╔═╡ 82292167-692a-442d-a6f5-9459fb1d013d
begin
data016 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_inf.α) + 1,  length(holstein_3d_alpha_0to12_beta_inf.ω) + 1)
data016[1, 1] = "κ(m₀ω₀²)|α↓|ω→"
data016[1, 2:end] = holstein_3d_alpha_0to12_beta_inf.ω
data016[2:end, 1] = 0.1:0.1:12
data016[2:end, 2:end] = reverse(ustrip(polaron_spring(holstein_3d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-spring-alpha-0to12-beta-inf.dat", data016)
data016
end

# ╔═╡ e0eddc80-9172-4d9f-af2e-3907d7af2b43
begin
data017 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_inf.α) + 1,  length(holstein_3d_alpha_0to12_beta_inf.ω) + 1)
data017[1, 1] = "R(r₀)|α↓|ω→"
data017[1, 2:end] = holstein_3d_alpha_0to12_beta_inf.ω
data017[2:end, 1] = 0.1:0.1:12
data017[2:end, 2:end] = reverse(ustrip(polaron_radius(holstein_3d_alpha_0to12_beta_inf)[:,:,1,1]),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-radius-alpha-0to12-beta-inf.dat", data017)
data017
end

# ╔═╡ 92b52583-b5db-43cd-8ee9-06098deece54
md"""
# 1D, 2D and 3D Holstein α = 0.1:0.1:12, ω = 1.0, T = LogRange(0.0625,32,128) ħω₀/kB
"""

# ╔═╡ 4589e7b4-5bfe-437f-927b-35626eb52b84
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 17121s
begin
holstein_1d_alpha_0to12_beta_003125to16 = holstein(12.0:-0.1:0.1, 1.0, 1.0, LogRange(16.0, 0.03125, 127), eps(); dims = 1)
save_holstein(holstein_1d_alpha_0to12_beta_003125to16, "data/holstein/variational/model/holstein-1d-alpha-0to12-beta-003125to16")
end
  ╠═╡ =#

# ╔═╡ ae307e8c-13c4-4d5b-8f99-4c7242df1565
holstein_1d_alpha_0to12_beta_003125to16 = load_holstein("data/holstein/variational/model/holstein-1d-alpha-0to12-beta-003125to16.jld")

# ╔═╡ 808e75f2-7956-4eef-bec8-0564c8750aec
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 5893s
begin
holstein_2d_alpha_0to12_beta_003125to16 = holstein(12.0:-0.1:0.1, 1.0, 1.0, LogRange(16.0, 0.03125, 127), eps(); dims = 2)
save_holstein(holstein_2d_alpha_0to12_beta_003125to16, "data/holstein/variational/model/holstein-2d-alpha-0to12-beta-003125to16")
end
  ╠═╡ =#

# ╔═╡ a9de638c-c8f7-438a-89fc-e362d7e725cb
holstein_2d_alpha_0to12_beta_003125to16 = load_holstein("data/holstein/variational/model/holstein-2d-alpha-0to12-beta-003125to16.jld")

# ╔═╡ 98274295-e983-4d0c-8c80-996c3d62ad37
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 2245s
begin
holstein_3d_alpha_0to12_beta_003125to16 = holstein(12.0:-0.1:0.1, 1.0, 1.0, LogRange(16.0, 0.03125, 127), eps(); dims = 3)
save_holstein(holstein_3d_alpha_0to12_beta_003125to16, "data/holstein/variational/model/holstein-3d-alpha-0to12-beta-003125to16")
end
  ╠═╡ =#

# ╔═╡ 966078d2-8692-4a6e-89b2-a227ee10cb0a
holstein_3d_alpha_0to12_beta_003125to16 = load_holstein("data/holstein/variational/model/holstein-3d-alpha-0to12-beta-003125to16.jld")

# ╔═╡ 8810d4e6-4ce8-4c29-af10-e34d7470a062
md"""
## Energy ħω₀
"""

# ╔═╡ a24c37d2-21d9-4493-a176-72bd81fb1b3a
md"""
### 1D
"""

# ╔═╡ 265d8a14-14a1-49a0-8beb-b22d423a5e9b
begin
data3 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_1d_alpha_0to12_beta_003125to16.α) + 1)
data3[1, 1] = "E(ħω₀)|T(ħω₀/kB)↓|α→"
data3[2:end, 1] = ustrip(1 ./ holstein_1d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data3[1, 2:end] = 0.1:0.1:12
data3[2:end, 2:end] = reverse(ustrip(holstein_1d_alpha_0to12_beta_003125to16.E[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat", data3)
data3
end

# ╔═╡ 6a881a63-e422-4389-b25b-23a2c6d82d29
md"""
### 2D
"""

# ╔═╡ c8d1b9fc-9c07-4ebb-8694-05f9586943eb
begin
data4 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data4[1, 1] = "E(ħω₀)|T(ħω₀/kB)↓|α→"
data4[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data4[1, 2:end] = 0.1:0.1:12
data4[2:end, 2:end] = reverse(ustrip(holstein_2d_alpha_0to12_beta_003125to16.E[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat", data4)
data4
end

# ╔═╡ f134f1d4-33c1-48c8-89ec-0f9f2d457e3e
md"""
### 3D
"""

# ╔═╡ 99faadff-2dd4-46ce-9a8c-d8faf850c9b7
begin
data5 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_3d_alpha_0to12_beta_003125to16.α) + 1)
data5[1, 1] = "E(ħω₀)|T(ħω₀/kB)↓|α→"
data5[2:end, 1] = ustrip(1 ./ holstein_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data5[1, 2:end] = 0.1:0.1:12
data5[2:end, 2:end] = reverse(ustrip(holstein_3d_alpha_0to12_beta_003125to16.E[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-3d-E-alpha-0to12-temp-00625to32.dat", data5)
data5
end

# ╔═╡ 958e1299-defb-410d-ac85-c176f621462c
md"""
## v ω₀
"""

# ╔═╡ 3bad3717-c2b5-4b2b-aa07-893c501c4215
md"""
### 1D
"""

# ╔═╡ 6ed07831-33d8-4d7f-acf9-41984568fb47
begin
data6 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_1d_alpha_0to12_beta_003125to16.α) + 1)
data6[1, 1] = "v(ω₀)|T(ħω₀/kB)↓|α→"
data6[2:end, 1] = ustrip(1 ./ holstein_1d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data6[1, 2:end] = 0.1:0.1:12
data6[2:end, 2:end] = reverse(ustrip(holstein_1d_alpha_0to12_beta_003125to16.v[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-1d-v-alpha-0to12-temp-00625to32.dat", data6)
data6
end

# ╔═╡ caa54305-0f13-47a4-850d-a3b14d975ad1
md"""
### 2D
"""

# ╔═╡ e5d9f83f-1c5d-4e0a-9842-554846c227dc
begin
data7 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data7[1, 1] = "v(ω₀)|T(ħω₀/kB)↓|α→"
data7[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data7[1, 2:end] = 0.1:0.1:12
data7[2:end, 2:end] = reverse(ustrip(holstein_2d_alpha_0to12_beta_003125to16.v[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat", data7)
data7
end

# ╔═╡ 13b93723-dcd3-4bc5-9935-30635ec2bbab
md"""
### 3D
"""

# ╔═╡ 1e306f70-be2b-4704-aa70-9d94d85d8859
begin
data8 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_3d_alpha_0to12_beta_003125to16.α) + 1)
data8[1, 1] = "v(ω₀)|T(ħω₀/kB)↓|α→"
data8[2:end, 1] = ustrip(1 ./ holstein_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data8[1, 2:end] = 0.1:0.1:12
data8[2:end, 2:end] = reverse(ustrip(holstein_3d_alpha_0to12_beta_003125to16.v[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-3d-v-alpha-0to12-temp-00625to32.dat", data8)
data8
end

# ╔═╡ c417194c-084c-440e-9bed-a5d7b9219ac8
md"""
## w ω₀
"""

# ╔═╡ e282d764-85fc-4f39-9d9e-2922353614b1
md"""
### 1D
"""

# ╔═╡ c1b0e649-e72a-487e-8475-349b9bfea8af
begin
data9 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_1d_alpha_0to12_beta_003125to16.α) + 1)
data9[1, 1] = "w(ω₀)|T(ħω₀/kB)↓|α→"
data9[2:end, 1] = ustrip(1 ./ holstein_1d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data9[1, 2:end] = 0.1:0.1:12
data9[2:end, 2:end] = reverse(ustrip(holstein_1d_alpha_0to12_beta_003125to16.w[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-1d-w-alpha-0to12-temp-00625to32.dat", data9)
data9
end

# ╔═╡ 2f330d5d-ebd5-4e56-b9d6-5e4dc8c48ac4
md"""
### 2D
"""

# ╔═╡ 670e50ac-afea-4e62-aef7-9bd12a8541d5
begin
data10 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data10[1, 1] = "w(ω₀)|T(ħω₀/kB)↓|α→"
data10[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data10[1, 2:end] = 0.1:0.1:12
data10[2:end, 2:end] = reverse(ustrip(holstein_2d_alpha_0to12_beta_003125to16.w[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-2d-w-alpha-0to12-temp-00625to32.dat", data10)
data10
end

# ╔═╡ e8c3adce-01b9-4ffe-a946-585cce2e4ae4
md"""
### 3D
"""

# ╔═╡ c3a3a96a-cc7c-4f13-9a01-14473007fcd9
begin
data11 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_3d_alpha_0to12_beta_003125to16.α) + 1)
data11[1, 1] = "w(ω₀)|T(ħω₀/kB)↓|α→"
data11[2:end, 1] = ustrip(1 ./ holstein_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data11[1, 2:end] = 0.1:0.1:12
data11[2:end, 2:end] = reverse(ustrip(holstein_3d_alpha_0to12_beta_003125to16.w[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-3d-w-alpha-0to12-temp-00625to32.dat", data11)
data11
end

# ╔═╡ 59947a9f-4599-4af8-a9c5-1f5f249ddc45
md"""
## M = v²/w² m₀
"""

# ╔═╡ 2127b843-5671-4f0f-afb4-0fa067398d36
md"""
### 1D
"""

# ╔═╡ 8288ed1a-d488-4f9b-aa3d-ed49c9b3844a
begin
data12 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_1d_alpha_0to12_beta_003125to16.α) + 1)
data12[1, 1] = "M(m₀)|T(ħω₀/kB)↓|α→"
data12[2:end, 1] = ustrip(1 ./ holstein_1d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data12[1, 2:end] = 0.1:0.1:12
data12[2:end, 2:end] = reverse(ustrip(polaron_mass(holstein_1d_alpha_0to12_beta_003125to16)[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-1d-mass-alpha-0to12-temp-00625to32.dat", data12)
data12
end

# ╔═╡ 1b6463e2-d744-4d7a-95fe-9d09ecd65e7f
md"""
### 2D
"""

# ╔═╡ 44c47eaa-bef8-4c50-bacd-b10866558caa
begin
data13 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data13[1, 1] = "M(m₀)|T(ħω₀/kB)↓|α→"
data13[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data13[1, 2:end] = 0.1:0.1:12
data13[2:end, 2:end] = reverse(ustrip(polaron_mass(holstein_2d_alpha_0to12_beta_003125to16)[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-2d-mass-alpha-0to12-temp-00625to32.dat", data13)
data13
end

# ╔═╡ 90a2c010-0650-4952-96d7-f8b03b7e8221
md"""
### 3D
"""

# ╔═╡ 9f1127a5-e401-4630-a429-71a862a54b44
begin
data14 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_3d_alpha_0to12_beta_003125to16.α) + 1)
data14[1, 1] = "M(m₀)|T(ħω₀/kB)↓|α→"
data14[2:end, 1] = ustrip(1 ./ holstein_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data14[1, 2:end] = 0.1:0.1:12
data14[2:end, 2:end] = reverse(ustrip(polaron_mass(holstein_3d_alpha_0to12_beta_003125to16)[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-3d-mass-alpha-0to12-temp-00625to32.dat", data14)
data14
end

# ╔═╡ 367c3815-0a85-493e-93bf-1d567ed46fe0
md"""
## κ = v²-w² m₀ω₀²
"""

# ╔═╡ c90ce2ea-38a9-4e1d-9804-07da5b9e02c9
md"""
### 1D
"""

# ╔═╡ 210ca750-f6ab-4664-8059-acfb1a4e7fa3
begin
data15 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_1d_alpha_0to12_beta_003125to16.α) + 1)
data15[1, 1] = "κ(m₀ω₀²)|T(ħω₀/kB)↓|α→"
data15[2:end, 1] = ustrip(1 ./ holstein_1d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data15[1, 2:end] = 0.1:0.1:12
data15[2:end, 2:end] = reverse(ustrip(polaron_spring(holstein_1d_alpha_0to12_beta_003125to16)[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat", data15)
data15
end

# ╔═╡ edb19919-7952-4578-93dd-78c1373d8036
md"""
### 2D
"""

# ╔═╡ a54fbca9-2752-45d1-9ded-06218ad89afe
begin
data16 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data16[1, 1] = "κ(m₀ω₀²)|T(ħω₀/kB)↓|α→"
data16[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data16[1, 2:end] = 0.1:0.1:12
data16[2:end, 2:end] = reverse(ustrip(polaron_spring(holstein_2d_alpha_0to12_beta_003125to16)[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-2d-spring-alpha-0to12-temp-00625to32.dat", data16)
data16
end

# ╔═╡ 52ce12a6-e895-4059-87ab-a1a26c6e05ea
md"""
### 3D
"""

# ╔═╡ fa667323-3183-41d0-8958-1ff5cefe7e9f
begin
data17 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data17[1, 1] = "κ(m₀ω₀²)|T(ħω₀/kB)↓|α→"
data17[2:end, 1] = ustrip(1 ./ holstein_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data17[1, 2:end] = 0.1:0.1:12
data17[2:end, 2:end] = reverse(ustrip(polaron_spring(holstein_3d_alpha_0to12_beta_003125to16)[:,1,1,:])',dims=2)
writedlm("data/holstein/variational/model/holstein-3d-spring-alpha-0to12-temp-00625to32.dat", data17)
data17
end

# ╔═╡ 76d0276b-47fd-4fb9-ae06-21c3d64a5f8a
md"""
## R = D/2 v/(v²-w²) coth(ħβv/2) r₀=√(ħ/m₀ω₀)
"""

# ╔═╡ 33c36b7b-73fa-4c0b-aa44-6487e64d7d71
md"""
### 1D
"""

# ╔═╡ cd07c2e0-6c9b-4737-bca9-98a0ed891178
begin
data18 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_1d_alpha_0to12_beta_003125to16.α) + 1)
data18[1, 1] = "R(r₀)|T(ħω₀/kB)↓|α→"
data18[2:end, 1] = ustrip(1 ./ holstein_1d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data18[1, 2:end] = 0.1:0.1:12
data18[2:end, 2:end] = reverse(ustrip(polaron_radius2(holstein_1d_alpha_0to12_beta_003125to16,dims=2)),dims=2)
writedlm("data/holstein/variational/model/holstein-1d-radius-alpha-0to12-temp-00625to32.dat", data18)
data18
end

# ╔═╡ 28e72ed4-aa97-482c-8dd4-35add756f650
md"""
### 2D
"""

# ╔═╡ e2b40339-f0da-4ee3-8017-1e18708b9175
begin
data19 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data19[1, 1] = "R(r₀)|T(ħω₀/kB)↓|α→"
data19[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data19[1, 2:end] = 0.1:0.1:12
data19[2:end, 2:end] = reverse(ustrip(polaron_radius2(holstein_2d_alpha_0to12_beta_003125to16,dims=2)),dims=2)
writedlm("data/holstein/variational/model/holstein-2d-radius-alpha-0to12-temp-00625to32.dat", data19)
data19
end

# ╔═╡ f5fb95ae-1be6-42cc-9c9f-8343590a1a4c
md"""
### 3D
"""

# ╔═╡ 658e98cb-5409-4b4b-9299-5d1ba8247cc1
begin
data20 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_3d_alpha_0to12_beta_003125to16.α) + 1)
data20[1, 1] = "R(r₀)|T(ħω₀/kB)↓|α→"
data20[2:end, 1] = ustrip(1 ./ holstein_3d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data20[1, 2:end] = 0.1:0.1:12
data20[2:end, 2:end] = reverse(ustrip(polaron_radius2(holstein_3d_alpha_0to12_beta_003125to16)),dims=2)
writedlm("data/holstein/variational/model/holstein-3d-radius-alpha-0to12-temp-00625to32.dat", data20)
data20
end

# ╔═╡ 790a9562-ce93-497b-9705-8fdb8985224c
md"""
## μ = 1/ℑΣ eω₀/m₀
"""

# ╔═╡ 1af79b7c-a9c1-4221-bcec-38db7af111af
md"""
### 1D
"""

# ╔═╡ 8bd1b064-00d7-495f-8f34-bfda55203a3f
begin
data21 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_1d_alpha_0to12_beta_003125to16.α) + 1)
data21[1, 1] = "μ(eω₀/m₀)|T(ħω₀/kB)↓|α→"
data21[2:end, 1] = ustrip(1 ./ holstein_1d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data21[1, 2:end] = 0.1:0.1:12
data21[2:end, 2:end] = reverse(ustrip(polaron_mobility(holstein_1d_alpha_0to12_beta_003125to16))',dims=2)
writedlm("data/holstein/variational/model/holstein-1d-mobility-alpha-0to12-temp-00625to32.dat", data21)
data21
end

# ╔═╡ d6aa74e2-d5d7-46d3-a7ca-27bba134f2d0
md"""
### 2D
"""

# ╔═╡ b3cca528-a7cb-442d-924e-d2fe83fad373
begin
data22 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_2d_alpha_0to12_beta_003125to16.α) + 1)
data22[1, 1] = "μ(eω₀/m₀)|T(ħω₀/kB)↓|α→"
data22[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data22[1, 2:end] = 0.1:0.1:12
data22[2:end, 2:end] = reverse(ustrip(polaron_mobility(holstein_2d_alpha_0to12_beta_003125to16))',dims=2)
writedlm("data/holstein/variational/model/holstein-2d-mobility-alpha-0to12-temp-00625to32.dat", data22)
data22
end

# ╔═╡ 3f1ad58b-71a0-4d65-8139-b79d37d765e0
md"""
### 3D
"""

# ╔═╡ 1f4a385a-3661-4bba-9d57-1615a0e4d60b
begin
data23 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_003125to16.β) + 1, length(holstein_3d_alpha_0to12_beta_003125to16.α) + 1)
data23[1, 1] = "μ(eω₀/m₀)|T(ħω₀/kB)↓|α→"
data23[2:end, 1] = ustrip(1 ./ holstein_2d_alpha_0to12_beta_003125to16.β ./ kB_pu)
data23[1, 2:end] = 0.1:0.1:12
data23[2:end, 2:end] = reverse(ustrip(polaron_mobility(holstein_3d_alpha_0to12_beta_003125to16))',dims=2)
writedlm("data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-temp-00625to32.dat", data23)
data23
end

# ╔═╡ 99354529-397a-41da-9bd0-c328a02025e1
md"""
# 1D, 2D and 3D Holstein α = 0.05:0.05:12, T = 0.48 ħω₀/kB, Ω = 0.01:0.01:30.0 ω₀
"""

# ╔═╡ a795f263-9236-4f8b-8e2f-21c8a0417738
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 19066s
begin
holstein_1d_alpha_0to12_beta_100_freq_0to30 = holstein(LinRange(12,0.05,240), 1.0, 1.0, 100.0, 0.01:0.01:30; dims = 1)
save_holstein(holstein_1d_alpha_0to12_beta_100_freq_0to30, "data/holstein/variational/model/holstein-1d-alpha-0to12-beta-100-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 9668b673-a459-491a-89f1-19979a9f0bac
holstein_1d_alpha_0to12_beta_100_freq_0to30 = load_holstein("data/holstein/variational/model/holstein-1d-alpha-0to12-beta-100-freq-0to30.jld")

# ╔═╡ 99859951-fce1-4224-9961-cb8bfb3a0ddb
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 2948s
begin
holstein_3d_alpha_0to12_beta_100_freq_0to30 = holstein(LinRange(12,0.05,240), 1.0, 1.0, 100.0, 0.01:0.01:30; dims = 3)
save_holstein(holstein_3d_alpha_0to12_beta_100_freq_0to30, "data/holstein/variational/model/holstein-3d-alpha-0to12-beta-100-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 9ae9e209-3db0-46e9-8e36-b4394a4d9417
holstein_3d_alpha_0to12_beta_100_freq_0to30 = load_holstein("data/holstein/variational/model/holstein-3d-alpha-0to12-beta-100-freq-0to30.jld")

# ╔═╡ 724fa1d6-cf9a-4c41-bff1-fc9108086d14
md"""
## Σ ω₀
"""

# ╔═╡ 5bf8f537-acac-4a5b-9fc3-f453d2ec2e29
md"""
### ℜΣ
"""

# ╔═╡ a8a79759-2277-45cf-9fd0-9b729cf7e465
begin
data24 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data24[1, 1] = "ℜΣ(ω₀)|α↓|Ω(ω₀)→"
data24[1, 2:end] = 0.01:0.01:30.0
data24[2:end, 1] = 0.05:0.05:12
data24[2:end, 2:end] = reverse(ustrip(real.(holstein_1d_alpha_0to12_beta_100_freq_0to30.Σ)),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-real-memory-alpha-0to12-beta-100-freq-0to30.dat", data24)
data24
end

# ╔═╡ c84cdaac-ceb0-4c16-884d-4d84df9d3b69
begin
data25 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data25[1, 1] = "ℜΣ(ω₀)|α↓|Ω(ω₀)→"
data25[1, 2:end] = 0.01:0.01:30.0
data25[2:end, 1] = 0.05:0.05:12
data25[2:end, 2:end] = reverse(ustrip(real.(holstein_2d_alpha_0to12_beta_100_freq_0to30.Σ)),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-real-memory-alpha-0to12-beta-100-freq-0to30.dat", data25)
data25
end

# ╔═╡ 6449a1f8-4faf-4216-80a1-9f2b63ca5d2e
begin
data26 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data26[1, 1] = "ℜΣ(ω₀)|α↓|Ω(ω₀)→"
data26[1, 2:end] = 0.01:0.01:30.0
data26[2:end, 1] = 0.05:0.05:12
data26[2:end, 2:end] = reverse(ustrip(real.(holstein_3d_alpha_0to12_beta_100_freq_0to30.Σ)),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-real-memory-alpha-0to12-beta-100-freq-0to30.dat", data26)
data26
end

# ╔═╡ d843c9af-7833-40ae-9391-90b3c55047ad
md"""
### ℑΣ
"""

# ╔═╡ 063a280d-e3d0-4db2-be8d-29486a846136
begin
data27 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data27[1, 1] = "ℑΣ(ω₀)|α↓|Ω(ω₀)→"
data27[1, 2:end] = 0.01:0.01:30.0
data27[2:end, 1] = 0.05:0.05:12
data27[2:end, 2:end] = reverse(ustrip(abs.(imag.(holstein_1d_alpha_0to12_beta_100_freq_0to30.Σ))),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat", data27)
data27
end

# ╔═╡ da93c141-d5e6-4af9-9a2b-ef6bb544abb3
begin
data28 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data28[1, 1] = "ℑΣ(ω₀)|α↓|Ω(ω₀)→"
data28[1, 2:end] = 0.01:0.01:30.0
data28[2:end, 1] = 0.05:0.05:12
data28[2:end, 2:end] = reverse(ustrip(abs.(imag.(holstein_2d_alpha_0to12_beta_100_freq_0to30.Σ))),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat", data28)
data28
end

# ╔═╡ f20bbb33-995c-46a7-871f-9090246a14cc
begin
data29 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data29[1, 1] = "ℑΣ(ω₀)|α↓|Ω(ω₀)→"
data29[1, 2:end] = 0.01:0.01:30.0
data29[2:end, 1] = 0.05:0.05:12
data29[2:end, 2:end] = reverse(ustrip(abs.(imag.(holstein_3d_alpha_0to12_beta_100_freq_0to30.Σ))),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat", data29)
data29
end

# ╔═╡ c7573414-d8b7-47db-8f01-9b3d7728117f
md"""
## σ = -im/(Ω+Σ) e²ω₀/m₀r₀²=e²ω₀²/ħ
"""

# ╔═╡ 9c8ecb33-095b-4526-b38f-f7ef7f3a50e9
md"""
### ℜσ
"""

# ╔═╡ bd7fa09a-039d-4ec1-8f50-9a52045462cf
begin
data30 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data30[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|α↓|Ω(ω₀)→"
data30[1, 2:end] = 0.01:0.01:30.0
data30[2:end, 1] = 0.05:0.05:12
data30[2:end, 2:end] = reverse(ustrip(abs.(real.(polaron_conductivity(holstein_1d_alpha_0to12_beta_100_freq_0to30)'))),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-real-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data30)
data30
end

# ╔═╡ a9f83dc8-c465-43cc-aef2-cfe82846bbfe
begin
data31 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data31[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|α↓|Ω(ω₀)→"
data31[1, 2:end] = 0.01:0.01:30.0
data31[2:end, 1] = 0.05:0.05:12
data31[2:end, 2:end] = reverse(ustrip(abs.(real.(polaron_conductivity(holstein_2d_alpha_0to12_beta_100_freq_0to30)'))),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-real-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data31)
data31
end

# ╔═╡ 2f816ee0-d6a2-4eb3-9027-a54461dd83bd
begin
data32 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data32[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|α↓|Ω(ω₀)→"
data32[1, 2:end] = 0.01:0.01:30.0
data32[2:end, 1] = 0.05:0.05:12
data32[2:end, 2:end] = reverse(ustrip(abs.(real.(polaron_conductivity(holstein_3d_alpha_0to12_beta_100_freq_0to30)'))),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-real-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data32)
data32
end

# ╔═╡ 962d8925-2b86-4af1-8971-21e9832c5faf
md"""
### ℑσ
"""

# ╔═╡ c31d6ebe-1d58-486c-a1aa-be9faf85b7e1
begin
data33 = Matrix{Any}(undef, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_1d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data33[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|α↓|Ω(ω₀)→"
data33[1, 2:end] = 0.01:0.01:30.0
data33[2:end, 1] = 0.05:0.05:12
data33[2:end, 2:end] = reverse(ustrip(imag.(polaron_conductivity(holstein_1d_alpha_0to12_beta_100_freq_0to30)')),dims=1)
writedlm("data/holstein/variational/model/holstein-1d-imag-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data33)
data33
end

# ╔═╡ 6882d7f1-0614-427d-9e01-f04741d3fa00
begin
data34 = Matrix{Any}(undef, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_2d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data34[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|α↓|Ω(ω₀)→"
data34[1, 2:end] = 0.01:0.01:30.0
data34[2:end, 1] = 0.05:0.05:12
data34[2:end, 2:end] = reverse(ustrip(imag.(polaron_conductivity(holstein_2d_alpha_0to12_beta_100_freq_0to30)')),dims=1)
writedlm("data/holstein/variational/model/holstein-2d-imag-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data34)
data34
end

# ╔═╡ e8a3b9af-24b0-4b6c-974f-5049ae1e4a43
begin
data35 = Matrix{Any}(undef, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.α) + 1, length(holstein_3d_alpha_0to12_beta_100_freq_0to30.Ω) + 1)
data35[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|α↓|Ω(ω₀)→"
data35[1, 2:end] = 0.01:0.01:30.0
data35[2:end, 1] = 0.05:0.05:12
data35[2:end, 2:end] = reverse(ustrip(imag.(polaron_conductivity(holstein_3d_alpha_0to12_beta_100_freq_0to30)')),dims=1)
writedlm("data/holstein/variational/model/holstein-3d-imag-conductivity-alpha-0to12-beta-100-freq-0to30.dat", data35)
data35
end

# ╔═╡ d83ae509-e1cd-4178-94bd-4037a68f2ac8
md"""
# 1D, 2D and 3D Holstein α = 1.0:7.0, T = LogRange(0.0625,32.0,128) ħω₀/kB, Ω = 0.01:0.01:30.0 ω₀
"""

# ╔═╡ e019f6ca-492c-4a22-a3e0-4ab46b2a0ca9
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 1973s 
for α in 1:4
	@eval $(Symbol("holstein_1d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(holstein(α, 1.0, 1.0, LogRange(0.0325,16.0,128), 0.01:0.01:30.0; dims = 1))
	save_holstein(eval(Symbol("holstein_1d_alpha_$(α)_beta_00325to16_freq_0to30")), "data/holstein/variational/model/holstein-1d-alpha-$(α)-beta-00325to16-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 6b9acd7e-661b-4580-9ce7-cc33526a5695
for α in 1:4
	@eval $(Symbol("holstein_1d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(load_holstein("data/holstein/variational/model/holstein-1d-alpha-$α-beta-00325to16-freq-0to30.jld"))
end

# ╔═╡ ed70cca0-be33-4fed-96d6-a21823514ccf
# ╠═╡ disabled = true
#=╠═╡
# thread 16 ~ 1915
for α in 1:5
	@eval $(Symbol("holstein_2d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(holstein(α, 1.0, 1.0, LogRange(0.0325,16.0,128), 0.01:0.01:30.0; dims = 2))
	save_holstein(eval(Symbol("holstein_2d_alpha_$(α)_beta_00325to16_freq_0to30")), "data/holstein/variational/model/holstein-2d-alpha-$(α)-beta-00325to16-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 1fd37ba6-1b71-4cf5-a69e-92184e7d7144
for α in 1:5
	@eval $(Symbol("holstein_2d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(load_holstein("data/holstein/variational/model/holstein-2d-alpha-$α-beta-00325to16-freq-0to30.jld"))
end

# ╔═╡ 4e4b75c3-582b-45b0-b389-49d9e2552817
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 2768s
for α in 1:7
	@eval $(Symbol("holstein_3d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(holstein(α, 1.0, 1.0, LogRange(0.0325,16.0,128), 0.01:0.01:30.0))
	save_holstein(eval(Symbol("holstein_3d_alpha_$(α)_beta_00325to16_freq_0to30")), "data/holstein/variational/model/holstein-3d-alpha-$(α)-beta-00325to16-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ 1d66f7a6-f8fa-478f-b007-956f24d047f7
for α in 1:7
	@eval $(Symbol("holstein_3d_alpha_$(α)_beta_00325to16_freq_0to30")) = $(load_holstein("data/holstein/variational/model/holstein-3d-alpha-$(α)-beta-00325to16-freq-0to30.jld"))
end

# ╔═╡ 9804abff-66b2-40af-8c16-a099e7600b73
md"""
## Σ ω₀
"""

# ╔═╡ 5fef1d73-a691-438b-8634-293c098f1522
md"""
### ℜΣ
"""

# ╔═╡ 100839bc-14ab-4c93-849b-b8ece778f189
for α in 1:4
	h = eval(Symbol("holstein_1d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data36 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data36[1, 1] = "ℜΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data36[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data36[1, 2:end] = 0.01:0.01:30
	data36[2:end, 2:end] = ustrip(real.(h.Σ))
	writedlm("data/holstein/variational/model/holstein-1d-real-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data36)
end

# ╔═╡ d06fbbf4-3f02-4b8f-a021-2e63c6ec4aef
for α in 1:5
	h = eval(Symbol("holstein_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data37 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data37[1, 1] = "ℜΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data37[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data37[1, 2:end] = 0.01:0.01:30
	data37[2:end, 2:end] = ustrip(real.(h.Σ))
	writedlm("data/holstein/variational/model/holstein-2d-real-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data37)
end

# ╔═╡ 9ce6fad2-9746-4605-b967-5447f8d0654d
for α in 1:7
	h = eval(Symbol("holstein_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data38 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data38[1, 1] = "ℜΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data38[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data38[1, 2:end] = 0.01:0.01:30
	data38[2:end, 2:end] = ustrip(real.(h.Σ))
	writedlm("data/holstein/variational/model/holstein-3d-real-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data38)
end

# ╔═╡ 45377ee4-b055-4ed4-8bb3-2b7321c4e709
md"""
### ℑΣ
"""

# ╔═╡ c9f99606-3824-4e7b-bb52-3d815c80a3f7
for α in 1:4
	h = eval(Symbol("holstein_1d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data39 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data39[1, 1] = "ℑΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data39[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data39[1, 2:end] = 0.01:0.01:30
	data39[2:end, 2:end] = ustrip(abs.(imag.(h.Σ)))
	writedlm("data/holstein/variational/model/holstein-1d-imag-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data39)
end

# ╔═╡ ca8ce8c3-d357-4913-868d-60944aee11e3
for α in 1:5
	h = eval(Symbol("holstein_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data40 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data40[1, 1] = "ℑΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data40[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data40[1, 2:end] = 0.01:0.01:30
	data40[2:end, 2:end] = ustrip(abs.(imag.(h.Σ)))
	writedlm("data/holstein/variational/model/holstein-2d-imag-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data40)
end

# ╔═╡ 20f9700f-3cb9-4197-8b9a-b9fe71cc748c
for α in 1:7
	h = eval(Symbol("holstein_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data41 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data41[1, 1] = "ℑΣ(ω₀)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data41[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data41[1, 2:end] = 0.01:0.01:30
	data41[2:end, 2:end] = ustrip(abs.(imag.(h.Σ)))
	writedlm("data/holstein/variational/model/holstein-3d-imag-memory-alpha-$α-temp-00625to32-freq-0to30.dat", data41)
end

# ╔═╡ b17f1e29-cc52-455b-a6ec-aec4dc70a372
md"""
## σ = -im/(Ω+Σ) e²ω₀/m₀r₀²=e²ω₀²/ħ
"""

# ╔═╡ 4d3645bb-cdf3-4e2c-9edd-c6c5be475107
md"""
### ℜσ
"""

# ╔═╡ 35ed15fa-4604-436e-9c8e-d4696f80afe4
for α in 1:4
	h = eval(Symbol("holstein_1d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data42 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data42[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data42[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data42[1, 2:end] = 0.01:0.01:30
	data42[2:end, 2:end] = ustrip(abs.(real.(polaron_conductivity(h))))'
	writedlm("data/holstein/variational/model/holstein-1d-real-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data42)
end

# ╔═╡ 8ede5bda-4291-41bd-b12c-7ec76df44ffc
for α in 1:5
	h = eval(Symbol("holstein_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data43 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data43[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data43[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data43[1, 2:end] = 0.01:0.01:30
	data43[2:end, 2:end] = ustrip(abs.(real.(polaron_conductivity(h))))'
	writedlm("data/holstein/variational/model/holstein-2d-real-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data43)
end

# ╔═╡ 84f94766-918f-40d6-a51a-2124f9e0df86
for α in 1:7
	h = eval(Symbol("holstein_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data44 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data44[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data44[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data44[1, 2:end] = 0.01:0.01:30
	data44[2:end, 2:end] = ustrip(abs.(real.(polaron_conductivity(h))))'
	writedlm("data/holstein/variational/model/holstein-3d-real-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data44)
end

# ╔═╡ fdc5be06-e1d0-4346-9402-0428fa451f8f
md"""
### ℑσ
"""

# ╔═╡ 36223675-c030-4c39-82a9-cab574387c9d
for α in 1:4
	h = eval(Symbol("holstein_1d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data45 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data45[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data45[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data45[1, 2:end] = 0.01:0.01:30
	data45[2:end, 2:end] = ustrip(imag.(polaron_conductivity(h)))'
	writedlm("data/holstein/variational/model/holstein-1d-imag-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data45)
end

# ╔═╡ a8b80879-0814-49ab-8f92-2255e1616af6
for α in 1:5
	h = eval(Symbol("holstein_2d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data46 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data46[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data46[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data46[1, 2:end] = 0.01:0.01:30
	data46[2:end, 2:end] = ustrip(imag.(polaron_conductivity(h)))'
	writedlm("data/holstein/variational/model/holstein-2d-imag-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data46)
end

# ╔═╡ 7f1fcea4-1933-442b-bdab-88c37e5b79b0
for α in 1:7
	h = eval(Symbol("holstein_3d_alpha_$(α)_beta_00325to16_freq_0to30"))
	data47 = Matrix{Any}(undef, length(h.β) + 1, length(h.Ω) + 1)
	data47[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|T(ħω₀/kB)↓|Ω(ω₀)→"
	data47[2:end, 1] = ustrip(1 ./ h.β ./ kB_pu)
	data47[1, 2:end] = 0.01:0.01:30
	data47[2:end, 2:end] = ustrip(imag.(polaron_conductivity(h)))'
	writedlm("data/holstein/variational/model/holstein-3d-imag-conductivity-alpha-$α-temp-00625to32-freq-0to30.dat", data47)
end

# ╔═╡ 0a53a5f4-30d5-425f-9691-6109ea9dbec3
holstein_2d_alpha_0to12_beta_100_freq_0to30 = load_holstein("data/holstein/variational/model/holstein-2d-alpha-0to12-beta-100-freq-0to30.jld")

# ╔═╡ 68457e83-c982-4326-94d2-2b5f9e60c6f7
# ╠═╡ disabled = true
#=╠═╡
# 16 threads ~ 3163s
begin
holstein_2d_alpha_0to12_beta_100_freq_0to30 = holstein(LinRange(12,0.05,240), 1.0, 1.0, 100.0, 0.01:0.01:30; dims = 2)
save_holstein(holstein_2d_alpha_0to12_beta_100_freq_0to30, "data/holstein/variational/model/holstein-2d-alpha-0to12-beta-100-freq-0to30")
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═bfd6e85c-48e2-11ef-2751-3147edce85d8
# ╠═35632de9-363d-459e-87b0-4bfaaccddb1f
# ╠═52e5925c-a82a-4799-97b4-a96a0798fa9b
# ╟─b12c5fa7-a588-4ce1-a4b5-1cf598380903
# ╟─4d74ec93-b41d-40af-b8e6-0652b2f10faf
# ╟─83427388-3d3f-4f05-867a-6984454cfdf9
# ╟─8564e569-76ed-439a-aa13-291538cfff13
# ╟─828f01c9-5726-4c8a-a5b1-bd6e83f29662
# ╟─6088c544-f44b-49c1-a503-ca881418e2f9
# ╟─8420455b-5772-413d-8e19-4e541e379142
# ╟─e03e239c-a77d-4591-aaa0-e02ce4cc0292
# ╟─6c2184de-0a93-43df-8b78-778b23844094
# ╟─4ee11b9d-f78c-437d-8e07-0c1fae31c2d5
# ╟─59777dfa-61f6-49fe-8b7b-5a68dd691c86
# ╟─057fde9d-dffc-4ca0-acde-813ba4c8966a
# ╟─45f46b37-e45d-4a01-ba51-9f00329876e6
# ╟─7b032f37-8a76-4db8-82c4-e7bf6719a031
# ╟─63b234ce-6b5f-4645-a3bb-04225f9575e7
# ╟─8cd17acc-9a62-45d6-8016-1af47643403a
# ╟─dd2c5cfd-b530-4380-9b5b-29739dc7217e
# ╟─9d95d5b9-1eee-4d7c-acd2-a2177f620fb7
# ╟─174468b5-1853-4eeb-858f-3a65f9c8a822
# ╟─c19ed360-aa81-4ed1-a5b6-0359b01e9a5f
# ╟─2cb035de-37aa-434b-b61e-a9106bac6b86
# ╟─7aec69dd-b635-4ae8-9be2-92f0fa8a6851
# ╟─be95e280-e5bd-4669-85f5-dbcb835b7b51
# ╟─82f79186-360b-4b78-a709-407bde56621c
# ╟─337a233b-3707-4103-b555-2e4fa94a1d14
# ╟─1f5dc90e-e61c-4213-8c02-09030b011a44
# ╟─90a6cf36-3b99-4841-bb6a-5c0192c11c37
# ╟─05caea00-de80-4357-83b3-8d884b8e0573
# ╟─9a91d18a-b49f-46b5-96e1-f577b142ceac
# ╟─477eaae2-731a-439a-9888-1d00f42943a8
# ╟─03d5c9cb-fe71-44c9-9838-f95283a3f874
# ╟─d04fe83b-a165-4c43-8316-906c1944fe5c
# ╟─0d60752c-f34a-4031-a1f4-57b8aa8c1717
# ╟─82292167-692a-442d-a6f5-9459fb1d013d
# ╟─e0eddc80-9172-4d9f-af2e-3907d7af2b43
# ╟─92b52583-b5db-43cd-8ee9-06098deece54
# ╟─4589e7b4-5bfe-437f-927b-35626eb52b84
# ╟─ae307e8c-13c4-4d5b-8f99-4c7242df1565
# ╟─808e75f2-7956-4eef-bec8-0564c8750aec
# ╟─a9de638c-c8f7-438a-89fc-e362d7e725cb
# ╟─98274295-e983-4d0c-8c80-996c3d62ad37
# ╟─966078d2-8692-4a6e-89b2-a227ee10cb0a
# ╟─8810d4e6-4ce8-4c29-af10-e34d7470a062
# ╟─a24c37d2-21d9-4493-a176-72bd81fb1b3a
# ╟─265d8a14-14a1-49a0-8beb-b22d423a5e9b
# ╟─6a881a63-e422-4389-b25b-23a2c6d82d29
# ╟─c8d1b9fc-9c07-4ebb-8694-05f9586943eb
# ╟─f134f1d4-33c1-48c8-89ec-0f9f2d457e3e
# ╟─99faadff-2dd4-46ce-9a8c-d8faf850c9b7
# ╟─958e1299-defb-410d-ac85-c176f621462c
# ╟─3bad3717-c2b5-4b2b-aa07-893c501c4215
# ╟─6ed07831-33d8-4d7f-acf9-41984568fb47
# ╟─caa54305-0f13-47a4-850d-a3b14d975ad1
# ╟─e5d9f83f-1c5d-4e0a-9842-554846c227dc
# ╟─13b93723-dcd3-4bc5-9935-30635ec2bbab
# ╟─1e306f70-be2b-4704-aa70-9d94d85d8859
# ╟─c417194c-084c-440e-9bed-a5d7b9219ac8
# ╟─e282d764-85fc-4f39-9d9e-2922353614b1
# ╟─c1b0e649-e72a-487e-8475-349b9bfea8af
# ╟─2f330d5d-ebd5-4e56-b9d6-5e4dc8c48ac4
# ╟─670e50ac-afea-4e62-aef7-9bd12a8541d5
# ╟─e8c3adce-01b9-4ffe-a946-585cce2e4ae4
# ╟─c3a3a96a-cc7c-4f13-9a01-14473007fcd9
# ╟─59947a9f-4599-4af8-a9c5-1f5f249ddc45
# ╟─2127b843-5671-4f0f-afb4-0fa067398d36
# ╟─8288ed1a-d488-4f9b-aa3d-ed49c9b3844a
# ╟─1b6463e2-d744-4d7a-95fe-9d09ecd65e7f
# ╟─44c47eaa-bef8-4c50-bacd-b10866558caa
# ╟─90a2c010-0650-4952-96d7-f8b03b7e8221
# ╟─9f1127a5-e401-4630-a429-71a862a54b44
# ╟─367c3815-0a85-493e-93bf-1d567ed46fe0
# ╟─c90ce2ea-38a9-4e1d-9804-07da5b9e02c9
# ╟─210ca750-f6ab-4664-8059-acfb1a4e7fa3
# ╟─edb19919-7952-4578-93dd-78c1373d8036
# ╟─a54fbca9-2752-45d1-9ded-06218ad89afe
# ╟─52ce12a6-e895-4059-87ab-a1a26c6e05ea
# ╟─fa667323-3183-41d0-8958-1ff5cefe7e9f
# ╟─76d0276b-47fd-4fb9-ae06-21c3d64a5f8a
# ╟─33c36b7b-73fa-4c0b-aa44-6487e64d7d71
# ╟─cd07c2e0-6c9b-4737-bca9-98a0ed891178
# ╟─28e72ed4-aa97-482c-8dd4-35add756f650
# ╟─e2b40339-f0da-4ee3-8017-1e18708b9175
# ╟─f5fb95ae-1be6-42cc-9c9f-8343590a1a4c
# ╟─658e98cb-5409-4b4b-9299-5d1ba8247cc1
# ╟─790a9562-ce93-497b-9705-8fdb8985224c
# ╟─1af79b7c-a9c1-4221-bcec-38db7af111af
# ╟─8bd1b064-00d7-495f-8f34-bfda55203a3f
# ╟─d6aa74e2-d5d7-46d3-a7ca-27bba134f2d0
# ╟─b3cca528-a7cb-442d-924e-d2fe83fad373
# ╟─3f1ad58b-71a0-4d65-8139-b79d37d765e0
# ╟─1f4a385a-3661-4bba-9d57-1615a0e4d60b
# ╟─99354529-397a-41da-9bd0-c328a02025e1
# ╟─a795f263-9236-4f8b-8e2f-21c8a0417738
# ╟─9668b673-a459-491a-89f1-19979a9f0bac
# ╟─68457e83-c982-4326-94d2-2b5f9e60c6f7
# ╟─0a53a5f4-30d5-425f-9691-6109ea9dbec3
# ╟─99859951-fce1-4224-9961-cb8bfb3a0ddb
# ╟─9ae9e209-3db0-46e9-8e36-b4394a4d9417
# ╟─724fa1d6-cf9a-4c41-bff1-fc9108086d14
# ╟─5bf8f537-acac-4a5b-9fc3-f453d2ec2e29
# ╟─a8a79759-2277-45cf-9fd0-9b729cf7e465
# ╟─c84cdaac-ceb0-4c16-884d-4d84df9d3b69
# ╟─6449a1f8-4faf-4216-80a1-9f2b63ca5d2e
# ╟─d843c9af-7833-40ae-9391-90b3c55047ad
# ╟─063a280d-e3d0-4db2-be8d-29486a846136
# ╟─da93c141-d5e6-4af9-9a2b-ef6bb544abb3
# ╟─f20bbb33-995c-46a7-871f-9090246a14cc
# ╟─c7573414-d8b7-47db-8f01-9b3d7728117f
# ╟─9c8ecb33-095b-4526-b38f-f7ef7f3a50e9
# ╟─bd7fa09a-039d-4ec1-8f50-9a52045462cf
# ╟─a9f83dc8-c465-43cc-aef2-cfe82846bbfe
# ╟─2f816ee0-d6a2-4eb3-9027-a54461dd83bd
# ╟─962d8925-2b86-4af1-8971-21e9832c5faf
# ╟─c31d6ebe-1d58-486c-a1aa-be9faf85b7e1
# ╟─6882d7f1-0614-427d-9e01-f04741d3fa00
# ╟─e8a3b9af-24b0-4b6c-974f-5049ae1e4a43
# ╟─d83ae509-e1cd-4178-94bd-4037a68f2ac8
# ╟─e019f6ca-492c-4a22-a3e0-4ab46b2a0ca9
# ╟─6b9acd7e-661b-4580-9ce7-cc33526a5695
# ╟─ed70cca0-be33-4fed-96d6-a21823514ccf
# ╟─1fd37ba6-1b71-4cf5-a69e-92184e7d7144
# ╟─4e4b75c3-582b-45b0-b389-49d9e2552817
# ╟─1d66f7a6-f8fa-478f-b007-956f24d047f7
# ╟─9804abff-66b2-40af-8c16-a099e7600b73
# ╟─5fef1d73-a691-438b-8634-293c098f1522
# ╠═100839bc-14ab-4c93-849b-b8ece778f189
# ╠═d06fbbf4-3f02-4b8f-a021-2e63c6ec4aef
# ╠═9ce6fad2-9746-4605-b967-5447f8d0654d
# ╟─45377ee4-b055-4ed4-8bb3-2b7321c4e709
# ╠═c9f99606-3824-4e7b-bb52-3d815c80a3f7
# ╠═ca8ce8c3-d357-4913-868d-60944aee11e3
# ╠═20f9700f-3cb9-4197-8b9a-b9fe71cc748c
# ╟─b17f1e29-cc52-455b-a6ec-aec4dc70a372
# ╟─4d3645bb-cdf3-4e2c-9edd-c6c5be475107
# ╠═35ed15fa-4604-436e-9c8e-d4696f80afe4
# ╠═8ede5bda-4291-41bd-b12c-7ec76df44ffc
# ╠═84f94766-918f-40d6-a51a-2124f9e0df86
# ╟─fdc5be06-e1d0-4346-9402-0428fa451f8f
# ╠═36223675-c030-4c39-82a9-cab574387c9d
# ╠═a8b80879-0814-49ab-8f92-2255e1616af6
# ╠═7f1fcea4-1933-442b-bdab-88c37e5b79b0
