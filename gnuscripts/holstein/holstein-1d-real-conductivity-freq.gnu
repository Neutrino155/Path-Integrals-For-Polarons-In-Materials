reset session

prefix = "../../plots/holstein/holstein-1d-real-conductivity-freq"
set key Left right top
set xlabel  "Frequency (ω₀)"
set ylabel  "Real Optical Absorption (ϵ₀cn)⁻¹"
set xrange [0:8]

plot    "../../data/holstein/variational/model/holstein-1d-real-conductivity-alpha-1to4-beta-100-freq-0to30.dat" u 1:2 w l t "{/Symbol a}=1", \
        "../../data/holstein/variational/model/holstein-1d-real-conductivity-alpha-1to4-beta-100-freq-0to30.dat" u 1:3 w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-1d-real-conductivity-alpha-1to4-beta-100-freq-0to30.dat" u 1:5 w l t "{/Symbol a}=3", \
        "../../data/holstein/variational/model/holstein-1d-real-conductivity-alpha-1to4-beta-100-freq-0to30.dat" u 1:6 w l t "{/Symbol a}=4"
     
load "../gnuplot-render.gpt"