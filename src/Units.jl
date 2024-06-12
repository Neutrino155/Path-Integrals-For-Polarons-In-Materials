# Units.jl

@unit   m0_pu  "m₀"    MassScale               1Unitful.me              false
@unit   e_pu   "e"     ElementaryCharge        1Unitful.q               false
@unit   ħ_pu   "ħ"     ReducedPlanckConstant   1Unitful.ħ               false
@unit   kB_pu  "kB"    BoltzmannConstant       1Unitful.k               false
@unit   ω0_pu  "ω₀"    PhononFrequency         1Unitful.THz2π           false
@unit   r0_pu  "r₀"    PolaronRadius           √(1ħ_pu/1m0_pu/1ω0_pu)   false
@unit   E0_pu  "ħω₀"   PhononEnergy            1ħ_pu*1ω0_pu             false

punit(x) = punit(dimension(x))

punit(x::Dimension{:Length})      = (r0_pu)^x.power
punit(x::Dimension{:Mass})        = (m0_pu)^x.power
punit(x::Dimension{:Time})        = (ħ_pu/E0_pu)^x.power
punit(x::Dimension{:Current})     = (e_pu*E0_pu/ħ_pu)^x.power
punit(x::Dimension{:Temperature}) = (E0_pu/kB_pu)^x.power

punit(::Dimension{D}) where D = throw(ArgumentError("No polaron unit defined for dimension $D."))

@generated punit(::Dimensions{N}) where N = prod(punit, N)
punit(::typeof(NoDims)) = NoUnits

for unit in (:(E0_pu), :(e_pu), :(ω0_pu), :(m0_pu), :(ħ_pu/r0_pu), :(ħ_pu), :(E0_pu/r0_pu),
    :(E0_pu/(e_pu*r0_pu)), :(ħ_pu/(e_pu*r0_pu^2)), :(E0_pu/e_pu), :(e_pu^2/(r0_pu*E0_pu)))
    @eval punit(::typeof(dimension($unit))) = $unit
end

puconvert(x) = uconvert(punit(x), x)
puconvert(u::Units, x::Number) = uconvert(u, x*punit(u))

pustrip(x) = ustrip(puconvert(x))
