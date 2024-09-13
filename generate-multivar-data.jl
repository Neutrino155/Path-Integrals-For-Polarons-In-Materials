### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 711042e8-6520-11ef-1758-491cd566c023
using Pkg, Unitful, Plots, DelimitedFiles

# ╔═╡ 07d58ada-9abe-47cf-a7fb-96b96b700b10
Pkg.activate(".")

# ╔═╡ 85d24807-1d4c-427b-9d4b-12d64f3464dc
using PolaronThesis

# ╔═╡ a6bae188-d450-45fa-a922-ea865628ba27
polaron_conductivity(frohlich) = -im ./ (frohlich.Ω .+ frohlich.Σ')

# ╔═╡ bc56571f-e612-4fc0-89fd-626f604cfd45
md"""
# Load Data
"""

# ╔═╡ 2027bdc8-95c0-4cfd-ad1c-c521f6d3e9f6
# ╠═╡ disabled = true
#=╠═╡
diagmc = readdlm("data/frohlich/diagmc/frohlich_mobility_freq_0to20_alpha_6_temp_05.txt")
  ╠═╡ =#

# ╔═╡ 4b6c74a7-adf4-4fdc-b3b8-0fcf724c7a97
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_inf_N_2 = load_frohlich("data/frohlich/variational/multivar/frohlich_3d_alpha_0to12_beta_inf_N_2.jld")
  ╠═╡ =#

# ╔═╡ e64f1c29-f756-4e9e-8885-3bde4478804b
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_inf_N_3 = load_frohlich("data/frohlich/variational/multivar/frohlich_3d_alpha_0to12_beta_inf_N_3.jld")
  ╠═╡ =#

# ╔═╡ 05ad9aee-416b-4504-afd4-83f7798dc4f9
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_inf_N_4 = load_frohlich("data/frohlich/variational/multivar/frohlich_3d_alpha_0to12_beta_inf_N_4.jld")
  ╠═╡ =#

# ╔═╡ b81a7e50-0645-4e0a-a539-974ea1947579
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_inf_N_5 = load_frohlich("data/frohlich/variational/multivar/frohlich_3d_alpha_0to12_beta_inf_N_5.jld")
  ╠═╡ =#

# ╔═╡ 88af9472-bc61-456c-a408-aea0524c4e0a
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_inf_N_6 = load_frohlich("data/frohlich/variational/multivar/frohlich_3d_alpha_0to12_beta_inf_N_6.jld")
  ╠═╡ =#

# ╔═╡ c23a3621-e8ce-4852-8caf-f30cebfa6960
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_003125to16_N_2 = load_frohlich("data/frohlich/variational/multivar/frohlich_3d_alpha_0to12_beta_003125to16_N_2.jld")
  ╠═╡ =#

# ╔═╡ 5406d718-81ab-42e1-93a7-41e0158217c6
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_003125to16_N_3 = load_frohlich("data/frohlich/variational/multivar/frohlich_3d_alpha_0to12_beta_003125to16_N_3.jld")
  ╠═╡ =#

# ╔═╡ d4609cdc-8457-4f5d-a254-b61d949b0b84
Ωrange = LinRange(0.01, 30, 3000)

# ╔═╡ 8838e520-b35c-4fdb-9e9e-3f8519b89785
# ╠═╡ disabled = true
#=╠═╡
# 19289s
begin
	frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2 = frohlich(frohlich_3d_alpha_0to12_beta_inf_N_2, Ωrange)
	save_frohlich(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2, "data/frohlich/variational/multivar/frohlich-3d-alpha-0to12-beta-inf-freq-0to30-N-2")
end
  ╠═╡ =#

# ╔═╡ 9f176c24-b300-4738-ac3e-1eecfe752d07
# ╠═╡ disabled = true
#=╠═╡
frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2 = load_frohlich("data/frohlich/variational/multivar/frohlich-3d-alpha-0to12-beta-inf-freq-0to30-N-2.jld")
  ╠═╡ =#

# ╔═╡ fcd11b51-124c-4513-86ae-3fd91afae836
# ╠═╡ disabled = true
#=╠═╡
begin
data1 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.α) + 1)
data1[1, 1] = "ℜΣ(ω₀)|Ω(ω₀)↓|α→"
data1[2:end, 1] = 0.01:0.01:30.0
data1[1, 2:end] = 0.1:0.1:12.0
data1[2:end, 2:end] = reverse(ustrip(real.(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.Σ)'),dims=2)
writedlm("data/frohlich/variational/multivar/frohlich-3d-real-memory-alpha-0to12-beta-100-freq-0to30-N-2.dat", data1)
data1
end
  ╠═╡ =#

# ╔═╡ e5067fa5-ac6e-4092-a5b5-77dc57fd0946
# ╠═╡ disabled = true
#=╠═╡
begin
data2 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.α) + 1)
data2[1, 1] = "ℑΣ(ω₀)|Ω(ω₀)↓|α→"
data2[2:end, 1] = 0.01:0.01:30.0
data2[1, 2:end] = 0.1:0.1:12.0
data2[2:end, 2:end] = reverse(abs.(ustrip(imag.(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.Σ))'),dims=2)
writedlm("data/frohlich/variational/multivar/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30-N-2.dat", data2)
data2
end
  ╠═╡ =#

# ╔═╡ c51ff0e7-8d5e-40fe-8c31-742995d1a267
# ╠═╡ disabled = true
#=╠═╡
begin
data3 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.α) + 1)
data3[1, 1] = "ℜσ(e²ω₀/m₀r₀²)|Ω(ω₀)↓|α→"
data3[2:end, 1] = 0.01:0.01:30.0
data3[1, 2:end] = 0.1:0.1:12.0
data3[2:end, 2:end] = reverse(ustrip(abs.(real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2)))),dims=2)
writedlm("data/frohlich/variational/multivar/frohlich-3d-real-conductivity-alpha-0to12-beta-100-freq-0to30-N-2.dat", data3)
data3
end
  ╠═╡ =#

# ╔═╡ a984c00a-fafa-438d-a793-0a6490e41863
# ╠═╡ disabled = true
#=╠═╡
begin
data4 = Matrix{Any}(undef, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.Ω) + 1, length(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2.α) + 1)
data4[1, 1] = "ℑσ(e²ω₀/m₀r₀²)|Ω(ω₀)↓|α→"
data4[2:end, 1] = 0.01:0.01:30.0
data4[1, 2:end] = 0.1:0.1:12.0
data4[2:end, 2:end] = reverse(ustrip(imag.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_inf_freq_0to30_N_2))),dims=2)
writedlm("data/frohlich/variational/multivar/frohlich-3d-imag-conductivity-alpha-0to12-beta-100-freq-0to30-N-2.dat", data4)
data4
end
  ╠═╡ =#

# ╔═╡ 2fcf096c-5ef7-4b7f-b7d3-70e33961f9cb
αrange = [0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]

# ╔═╡ 44e69331-2672-46a1-968d-2ec1435306e6
frohlich_3d_alpha_0to12_beta_100_N_1 = frohlich(αrange, 1.0, 100.0)

# ╔═╡ 7d8b473f-0290-4724-8deb-a390d72ce74d
frohlich_3d_alpha_0to12_beta_100_N_1_0to30 = frohlich(frohlich_3d_alpha_0to12_beta_100_N_1, 0.1:0.1:10.0)

# ╔═╡ 78d660b5-eb45-41a7-b292-4432e3aaae7e
plot(0.1:0.1:10.0, real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_N_1_0to30)), ylims=(0,0.5))

# ╔═╡ 1554c91f-d325-4800-bc17-020491ed5c2e
frohlich_3d_alpha_0to12_beta_100_N_2 = frohlich(αrange, 1.0, 100.0; v_guesses = [6, 4], w_guesses = [5, 3], upper = 100.0)

# ╔═╡ 271f0b4e-224b-45f8-bd8b-5e8b9ae616be
frohlich_3d_alpha_0to12_beta_100_N_2_0to30 = frohlich(frohlich_3d_alpha_0to12_beta_100_N_2, 0.1:0.1:10.0)

# ╔═╡ ef3eaac1-f9d7-4e8a-bd71-11c6992985f6
plot(0.1:0.1:10.0, real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_N_2_0to30)), ylims=(0,0.5))

# ╔═╡ de5cf838-1f44-461e-b7c6-9c31eb7f1c40
frohlich_3d_alpha_0to12_beta_100_N_3 = frohlich(αrange, 1.0, 100.0; v_guesses = [8,6, 4], w_guesses = [7,5, 3], upper = 200.0)

# ╔═╡ 3bbac0ae-c86c-4e0a-a460-d63b22991a44
frohlich_3d_alpha_0to12_beta_100_N_3_0to30 = frohlich(frohlich_3d_alpha_0to12_beta_100_N_3, 0.1:0.1:10.0)

# ╔═╡ 14c47316-b5fd-4d6c-a7d6-6ebda865f989
plot(0.1:0.1:10.0, real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_N_3_0to30)), ylims=(0,0.5))

# ╔═╡ ad676377-f6be-46a3-81fd-293f0085d2c7
frohlich_3d_alpha_0to12_beta_100_N_4 = frohlich(αrange, 1.0, 100.0; v_guesses = [10,8,6,4], w_guesses = [9,7,5,3], upper = 500.0)

# ╔═╡ b10fb360-a44f-4372-bd7a-5c6a794a5be7
frohlich_3d_alpha_0to12_beta_100_N_4_0to30 = frohlich(frohlich_3d_alpha_0to12_beta_100_N_4, 0.1:0.1:10.0)

# ╔═╡ 03ac42ff-7a2c-4b61-a6d3-4500c4913080
plot(0.1:0.1:10.0, real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_N_4_0to30)), ylims=(0,0.5))

# ╔═╡ 8bf786ef-4d7a-4a8f-bff0-245b3e31a30e
frohlich_3d_alpha_0to12_beta_100_N_5 = frohlich(αrange, 1.0, 100.0; v_guesses = [12,10,8,6,4], w_guesses = [11,9,7,5,3], upper = 1000.0)

# ╔═╡ 88e8b2c6-0a36-4ce8-a36b-03a938209fb6
frohlich_3d_alpha_0to12_beta_100_N_5_0to30 = frohlich(frohlich_3d_alpha_0to12_beta_100_N_5, 0.1:0.1:10.0)

# ╔═╡ 404240e9-5f66-445e-a0ff-bfd6568fd23f
plot(0.1:0.1:10.0, real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_N_5_0to30)), ylims=(0,0.5))

# ╔═╡ 9518e064-89f1-4205-b88c-c7894c48d4b9
frohlich_3d_alpha_0to12_beta_100_N_6 = frohlich(αrange, 1.0, 100.0; v_guesses = [14,12,10,8,6,4], w_guesses = [13,11,9,7,5,3], upper = 2000.0)

# ╔═╡ 56c43383-0a6d-4bae-bdee-a3684848d169
frohlich_3d_alpha_0to12_beta_100_N_6_0to30 = frohlich(frohlich_3d_alpha_0to12_beta_100_N_6, 0.1:0.1:10.0)

# ╔═╡ 56304888-bd8d-43a4-996f-7f2967c06936
plot(0.1:0.1:10.0, real.(polaron_conductivity(frohlich_3d_alpha_0to12_beta_100_N_6_0to30)), ylims=(0,0.5))

# ╔═╡ 4c5b383f-1db5-40ae-8ab0-e3f7fbff501c
frohlich_3d_alpha_0to12_beta_100_N_7 = frohlich(αrange, 1.0, 100.0; v_guesses = [16,14,12,10,8,6,4], w_guesses = [15,13,11,9,7,5,3], upper = 2000.0)

# ╔═╡ d37d9cf6-037d-413a-961b-83bde3642892
begin
	plot(αrange, reduce_array(frohlich_3d_alpha_0to12_beta_100_N_1.E .- frohlich_3d_alpha_0to12_beta_100_N_2.E), yaxis=:log)
	plot!(αrange, reduce_array(frohlich_3d_alpha_0to12_beta_100_N_2.E .- frohlich_3d_alpha_0to12_beta_100_N_3.E))
	plot!(αrange, reduce_array(frohlich_3d_alpha_0to12_beta_100_N_3.E .- frohlich_3d_alpha_0to12_beta_100_N_4.E))
	plot!(αrange, reduce_array(frohlich_3d_alpha_0to12_beta_100_N_4.E .- frohlich_3d_alpha_0to12_beta_100_N_5.E))
	plot!(αrange, reduce_array(frohlich_3d_alpha_0to12_beta_100_N_5.E .- frohlich_3d_alpha_0to12_beta_100_N_6.E))
	plot!(αrange, reduce_array(frohlich_3d_alpha_0to12_beta_100_N_6.E .- frohlich_3d_alpha_0to12_beta_100_N_7.E))
end

# ╔═╡ cfc13eed-2322-4052-8075-e96354c56287
reduce_array(frohlich_3d_alpha_0to12_beta_100_N_6.E .- frohlich_3d_alpha_0to12_beta_100_N_7.E)

# ╔═╡ Cell order:
# ╠═711042e8-6520-11ef-1758-491cd566c023
# ╟─07d58ada-9abe-47cf-a7fb-96b96b700b10
# ╠═85d24807-1d4c-427b-9d4b-12d64f3464dc
# ╠═a6bae188-d450-45fa-a922-ea865628ba27
# ╟─bc56571f-e612-4fc0-89fd-626f604cfd45
# ╠═2027bdc8-95c0-4cfd-ad1c-c521f6d3e9f6
# ╠═4b6c74a7-adf4-4fdc-b3b8-0fcf724c7a97
# ╠═e64f1c29-f756-4e9e-8885-3bde4478804b
# ╠═05ad9aee-416b-4504-afd4-83f7798dc4f9
# ╠═b81a7e50-0645-4e0a-a539-974ea1947579
# ╠═88af9472-bc61-456c-a408-aea0524c4e0a
# ╠═c23a3621-e8ce-4852-8caf-f30cebfa6960
# ╠═5406d718-81ab-42e1-93a7-41e0158217c6
# ╟─d4609cdc-8457-4f5d-a254-b61d949b0b84
# ╠═8838e520-b35c-4fdb-9e9e-3f8519b89785
# ╠═9f176c24-b300-4738-ac3e-1eecfe752d07
# ╠═fcd11b51-124c-4513-86ae-3fd91afae836
# ╠═e5067fa5-ac6e-4092-a5b5-77dc57fd0946
# ╠═c51ff0e7-8d5e-40fe-8c31-742995d1a267
# ╠═a984c00a-fafa-438d-a793-0a6490e41863
# ╠═2fcf096c-5ef7-4b7f-b7d3-70e33961f9cb
# ╠═44e69331-2672-46a1-968d-2ec1435306e6
# ╠═7d8b473f-0290-4724-8deb-a390d72ce74d
# ╠═78d660b5-eb45-41a7-b292-4432e3aaae7e
# ╠═1554c91f-d325-4800-bc17-020491ed5c2e
# ╠═271f0b4e-224b-45f8-bd8b-5e8b9ae616be
# ╠═ef3eaac1-f9d7-4e8a-bd71-11c6992985f6
# ╠═de5cf838-1f44-461e-b7c6-9c31eb7f1c40
# ╠═3bbac0ae-c86c-4e0a-a460-d63b22991a44
# ╠═14c47316-b5fd-4d6c-a7d6-6ebda865f989
# ╠═ad676377-f6be-46a3-81fd-293f0085d2c7
# ╠═b10fb360-a44f-4372-bd7a-5c6a794a5be7
# ╠═03ac42ff-7a2c-4b61-a6d3-4500c4913080
# ╠═8bf786ef-4d7a-4a8f-bff0-245b3e31a30e
# ╠═88e8b2c6-0a36-4ce8-a36b-03a938209fb6
# ╠═404240e9-5f66-445e-a0ff-bfd6568fd23f
# ╠═9518e064-89f1-4205-b88c-c7894c48d4b9
# ╠═56c43383-0a6d-4bae-bdee-a3684848d169
# ╠═56304888-bd8d-43a4-996f-7f2967c06936
# ╠═4c5b383f-1db5-40ae-8ab0-e3f7fbff501c
# ╠═d37d9cf6-037d-413a-961b-83bde3642892
# ╠═cfc13eed-2322-4052-8075-e96354c56287
