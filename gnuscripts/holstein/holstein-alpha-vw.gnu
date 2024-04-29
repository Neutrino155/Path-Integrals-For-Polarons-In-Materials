reset session

prefix = "../../plots/holstein/holstein-alpha-vw"
set key left
set grid
set xlabel  "α"
set ylabel  "v, w"
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/holstein/variational/model/holstein-1d-v-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "v 1d" lw 2, \
        "../../data/holstein/variational/model/holstein-1d-w-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "w 1d" lw 2, \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "v 2d" lw 2, \
        "../../data/holstein/variational/model/holstein-2d-w-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "w 2d" lw 2, \
        "../../data/holstein/variational/model/holstein-3d-v-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "v 3d" lw 2, \
        "../../data/holstein/variational/model/holstein-3d-w-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "w 3d" lw 2

load "../gnuplot-render.gpt"