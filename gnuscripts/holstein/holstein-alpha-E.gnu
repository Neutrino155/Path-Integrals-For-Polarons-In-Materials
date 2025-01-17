reset session

prefix = "../../plots/holstein/holstein-alpha-energy"
set key right
set grid
set xlabel  "α"
set ylabel  "Energy (ħω₀)"
set xrange [0:12]
set yrange [-18:0]
set xtics 0,1,12

plot    "../../data/frohlich/variational/model/frohlich-3d-alpha-0to12-beta-inf.dat" u 1:4 w l t "Fro", \
        "../../data/holstein/variational/model/holstein-1d-energy-alpha-0to12-beta-inf.dat" u 1:11 w l t "1D", \
        "../../data/holstein/variational/model/holstein-2d-energy-alpha-0to12-beta-inf.dat" u 1:11 w l t "2D", \
        "../../data/holstein/variational/model/holstein-3d-energy-alpha-0to12-beta-inf.dat" u 1:11 w l t "3D"
       
load "../gnuplot-render.gpt"