reset session

prefix = "../../plots/frohlich/frohlich-2d-alpha-0to12"
set key Left left top
set grid
set xlabel  "α"
set ylabel  "Frohlich Properties"
set xrange [0:12]
set yrange [4**-3:4**5]
set for [i=-3:5] ytics (sprintf("4^{%d}", i) 4**i)
set xtics 0,1,12
set logscale y 2

plot    "../../data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.dat" u 1:2 w l t "v", \
        "../../data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.dat" u 1:3 w l t "w", \
        "../../data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.dat" u 1:(-$4) w l t "-E", \
        "../../data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.dat" u 1:5 w l t "M", \
        "../../data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.dat" u 1:6 w l t "R", \
        "../../data/frohlich/variational/model/frohlich-2d-alpha-0to12-beta-inf.dat" u 1:7 w l t "κ"

load "../gnuplot-render.gpt"
