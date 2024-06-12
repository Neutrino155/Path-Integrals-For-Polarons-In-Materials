reset session

prefix = "../../plots/holstein/holstein-2d-imag-conductivity-freq"
set key Left right top
set xlabel  "Frequency (ω₀)"
set ylabel  "Holstein Imag Optical Absorption (ϵ₀cn)⁻¹"
set yrange [0:0.4]
set xrange [0:16]

plot    "../../data/holstein/variational/model/holstein-2d-imag-conductivity-alpha-1to5-beta-100-freq-0to30.dat" u 1:2 w l t "{/Symbol a}=1", \
        "../../data/holstein/variational/model/holstein-2d-imag-conductivity-alpha-1to5-beta-100-freq-0to30.dat" u 1:3 w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-2d-imag-conductivity-alpha-1to5-beta-100-freq-0to30.dat" u 1:5 w l t "{/Symbol a}=3", \
        "../../data/holstein/variational/model/holstein-2d-imag-conductivity-alpha-1to5-beta-100-freq-0to30.dat" u 1:6 w l t "{/Symbol a}=4", \
        "../../data/holstein/variational/model/holstein-2d-imag-conductivity-alpha-1to5-beta-100-freq-0to30.dat" u 1:7 w l t "{/Symbol a}=5"

load "../gnuplot-render.gpt"