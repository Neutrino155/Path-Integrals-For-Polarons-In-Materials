reset session

prefix = "../../plots/rubrene/rubrene-mobility-temp-1to400K"
set key right
set grid
set xlabel  "T (K)"
set ylabel  "Mobility μ (cm²V⁻¹s⁻¹)"
set xrange [1:400]
set xtics (1,50,100,150,200,250,300,350,400)
set logscale y 10
set yrange [10:1e4]

plot    "../../data/holstein/variational/Rubrene/holstein-rubrene-temp-1to400K.dat" u 1:7 w l t "H", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-temp-1to400K.dat" u 1:7 w l t "F"
       
load "../gnuplot-render.gpt"