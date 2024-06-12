reset session

prefix = "../../plots/rubrene/rubrene-vw-temp-1to400K"
set key bottom right
set grid
set xlabel  "T (K)"
set ylabel  "v, w (ω₀)"
set xrange [1:400]
set xtics (1,50,100,150,200,250,300,350,400)

plot    "../../data/holstein/variational/Rubrene/holstein-rubrene-temp-1to400K.dat" u 1:2 w l t "v_{H}", \
        "../../data/holstein/variational/Rubrene/holstein-rubrene-temp-1to400K.dat" u 1:3 w l t "w_{H}", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-temp-1to400K.dat" u 1:2 w l t "v_{F}", \
        "../../data/frohlich/variational/Rubrene/frohlich-rubrene-temp-1to400K.dat" u 1:3 w l t "w_{F}"
       
load "../gnuplot-render.gpt"