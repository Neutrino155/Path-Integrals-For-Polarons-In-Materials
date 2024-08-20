reset session

prefix = "../../plots/frohlich/frohlich-3d-imag-memory-freq"
set key right bottom
set grid
set xlabel  "Frequency (ω₀)"
set ylabel  "Frohlich Imag Memory (ω₀m₀⁻¹)"
set yrange [0.01:100]
set logscale y

plot    "../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($21 * $1) w l t "{/Symbol a}=1", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($41 * $1) w l t "{/Symbol a}=2", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($61 * $1) w l t "{/Symbol a}=3", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($81 * $1) w l t "{/Symbol a}=4", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($101 * $1) w l t "{/Symbol a}=5", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($121 * $1) w l t "{/Symbol a}=6", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:($141 * $1) w l t "{/Symbol a}=7"

load "../gnuplot-render.gpt"