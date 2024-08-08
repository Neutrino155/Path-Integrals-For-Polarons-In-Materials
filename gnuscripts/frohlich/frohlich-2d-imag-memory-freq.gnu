reset session

prefix = "../../plots/frohlich/frohlich-2d-imag-memory-freq"
set key Left right top
set xlabel  "Frequency (ω₀)"
set ylabel  "Frohlich Imag Memory (ω₀m₀⁻¹)"
set yrange [0:*]

plot    "../../data/frohlich/variational/model/frohlich-2d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:21 w l t "{/Symbol a}=1", \
        "../../data/frohlich/variational/model/frohlich-2d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:41 w l t "{/Symbol a}=2", \
        "../../data/frohlich/variational/model/frohlich-2d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:61 w l t "{/Symbol a}=3", \
        "../../data/frohlich/variational/model/frohlich-2d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:81 w l t "{/Symbol a}=4", \
        "../../data/frohlich/variational/model/frohlich-2d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat" u 1:101 w l t "{/Symbol a}=5"

load "../gnuplot-render.gpt"