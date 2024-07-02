reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-alpha-0to12-beta-inf-N-1to6"
set key Left left top
set grid
set xlabel  "α"
set ylabel  "Percentage Difference with N=1"
set xrange [0:12]
set yrange [0:0.015]
set xtics 0,1,12

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-0to12-beta-inf-N-1to6.dat" u 1:($7 - $2) w l t "N=6", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-0to12-beta-inf-N-1to6.dat" u 1:($2 - $6) w l t "N=5", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-0to12-beta-inf-N-1to6.dat" u 1:($2 - $5) w l t "N=4", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-0to12-beta-inf-N-1to6.dat" u 1:($2 - $4) w l t "N=3", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-0to12-beta-inf-N-1to6.dat" u 1:($2 - $3) w l t "N=2", \
        "../../data/frohlich/variational/multivar/dries-sels-alpha-0to15-beta-inf-N-inf.dat" u 1:2 pt 7 dt 1 t "Best"

load "../gnuplot-render.gpt"
