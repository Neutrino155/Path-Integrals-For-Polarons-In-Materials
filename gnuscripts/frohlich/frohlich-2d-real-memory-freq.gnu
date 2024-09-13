reset session

prefix = "../../plots/frohlich/frohlich-2d-real-memory-freq"
set key right top
set grid
set xlabel  "Frequency (ω₀)"
set ylabel  "Frohlich Real Memory (ω₀m₀⁻¹)"
set yrange [-25:25]

plot    "../../data/frohlich/variational/model/frohlich-2d-real-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($21) w l t "{/Symbol a}=1", \
        "../../data/frohlich/variational/model/frohlich-2d-real-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($41) w l t "{/Symbol a}=2", \
        "../../data/frohlich/variational/model/frohlich-2d-real-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($61) w l t "{/Symbol a}=3", \
        "../../data/frohlich/variational/model/frohlich-2d-real-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($81) w l t "{/Symbol a}=4", \
        "../../data/frohlich/variational/model/frohlich-2d-real-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($101) w l t "{/Symbol a}=5"

load "../gnuplot-render.gpt"