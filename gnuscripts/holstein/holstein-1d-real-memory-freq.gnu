reset session

prefix = "../../plots/holstein/holstein-1d-real-memory-freq"
set key Left right top
set xlabel  "Frequency (ω₀)"
set ylabel  "Holstein Real Memory (m/ω₀)"

plot    "../../data/holstein/variational/model/holstein-1d-real-memory-alpha-1to4-beta-100-freq-0to30.dat" u 1:2 w l t "{/Symbol a}=1", \
        "../../data/holstein/variational/model/holstein-1d-real-memory-alpha-1to4-beta-100-freq-0to30.dat" u 1:3 w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-1d-real-memory-alpha-1to4-beta-100-freq-0to30.dat" u 1:5 w l t "{/Symbol a}=3", \
        "../../data/holstein/variational/model/holstein-1d-real-memory-alpha-1to4-beta-100-freq-0to30.dat" u 1:6 w l t "{/Symbol a}=4"
     
load "../gnuplot-render.gpt"