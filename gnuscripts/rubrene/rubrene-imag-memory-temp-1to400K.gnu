reset session

prefix = "../../plots/rubrene/rubrene-imag-memory-temp-1to400K"
set key Left left
set grid
set xlabel  "T (K)"
set ylabel  "Imaginary Memory Imχ (m₀/ω₀)"
set xrange [1:400]
set xtics (1,50,100,150,200,250,300,350,400)
set logscale y 10
set yrange [1e-3:*]

plot    "../../data/holstein/variational/Rubrene/holstein-rubrene-temp-1to400K.dat" u 1:6 w l t "H", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-temp-1to400K.dat" u 1:6 w l t "F"
       
load "../gnuplot-render.gpt"