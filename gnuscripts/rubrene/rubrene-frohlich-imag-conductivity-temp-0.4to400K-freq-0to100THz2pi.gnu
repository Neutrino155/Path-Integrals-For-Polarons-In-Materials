reset session

prefix = "../../plots/rubrene/rubrene-frohlich-imag-conductivity-temp-0.4to400K-freq-0to100THz2pi"
set key top right
set grid
set xlabel  "Frequency (THz2π)"
set ylabel  "Imag Conductivity σ (μS)"
set xrange [0:100]
set xtics 0,10,100
set mxtics 2
set yrange [0:1000]

plot    "../../data/frohlich/variational/Rubrene/frohlich-rubrene-imag-conductivity-temp-0.4to400K-freq-0to100THz2pi.dat" u 1:(-$2) w l t "0.48 K", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-imag-conductivity-temp-0.4to400K-freq-0to100THz2pi.dat" u 1:(-$4) w l t "100 K", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-imag-conductivity-temp-0.4to400K-freq-0to100THz2pi.dat" u 1:(-$6) w l t "200 K", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-imag-conductivity-temp-0.4to400K-freq-0to100THz2pi.dat" u 1:(-$8) w l t "300 K", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-imag-conductivity-temp-0.4to400K-freq-0to100THz2pi.dat" u 1:(-$10) w l t "400 K"
       
load "../gnuplot-render.gpt"