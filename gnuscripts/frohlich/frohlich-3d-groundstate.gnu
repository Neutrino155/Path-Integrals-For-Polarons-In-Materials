reset session

prefix = "../../plots/frohlich/frohlich-3d-alpha-0to12"
set key Left left top
set grid
set xlabel  "α"
set ylabel  "Frohlich Properties" offset 4,0
set xrange [0:12]
set yrange [*:*]
set xtics 0,1,12
set lmargin 8
set logscale y 2
set for [i=-2:4] ytics add (sprintf("4^{%d}", i) 4**i)

plot    "../../data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat" u 1:2 w l t "v", \
        "../../data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat" u 1:3 w l t "w", \
        "../../data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat" u 1:(-$4) w l t "-E", \
        "../../data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat" u 1:5 w l t "M", \
        "../../data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat" u 1:6 w l t "R", \
        "../../data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat" u 1:7 w l t "κ"

load "../gnuplot-render.gpt"
