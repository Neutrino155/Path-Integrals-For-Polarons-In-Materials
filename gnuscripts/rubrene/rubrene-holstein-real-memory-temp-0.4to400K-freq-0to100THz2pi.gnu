reset session

prefix = "../../plots/rubrene/rubrene-holstein-real-memory-temp-0.4to400K-freq-0to10omega"
set key top left
set grid
set xlabel  "Frequency (THz2π)"
set ylabel  "Real Memory Reχ (ω₀/m₀)"
set mxtics 2
set xrange [0:50]

plot    "../../data/holstein/variational/Rubrene/holstein-rubrene-real-memory-temp-0.4to400K-freq-0to10omega.dat" u 1:2 w l t "0.48 K", \
        "../../data/holstein/variational/Rubrene/holstein-rubrene-real-memory-temp-0.4to400K-freq-0to10omega.dat" u 1:4 w l t "100 K", \
        "../../data/holstein/variational/Rubrene/holstein-rubrene-real-memory-temp-0.4to400K-freq-0to10omega.dat" u 1:6 w l t "200 K", \
        "../../data/holstein/variational/Rubrene/holstein-rubrene-real-memory-temp-0.4to400K-freq-0to10omega.dat" u 1:8 w l t "300 K", \
        "../../data/holstein/variational/Rubrene/holstein-rubrene-real-memory-temp-0.4to400K-freq-0to10omega.dat" u 1:10 w l t "400 K"
       
load "../gnuplot-render.gpt"