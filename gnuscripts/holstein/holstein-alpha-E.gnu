reset session

prefix = "../../plots/holstein/holstein-alpha-energy"
set key right
set grid
set xlabel  "α"
set ylabel  "Energy (J)"
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "1d", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "2d", \
        "../../data/holstein/variational/model/holstein-3d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "3d"
       
load "../gnuplot-render.gpt"