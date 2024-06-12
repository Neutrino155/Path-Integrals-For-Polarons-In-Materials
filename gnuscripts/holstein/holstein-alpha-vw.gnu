reset session

prefix = "../../plots/holstein/holstein-alpha-vw"
set key Left left
set grid
set xlabel  "α"
set ylabel  "v, w"
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/frohlich/variational/model/frohlich-alpha-0to12-beta-inf.dat" u 1:2 w l t "v_{F}", \
        "../../data/frohlich/variational/model/frohlich-alpha-0to12-beta-inf.dat" u 1:3 w l t "w_{F}", \
        "../../data/holstein/variational/model/holstein-1d-v-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "v_{1D}", \
        "../../data/holstein/variational/model/holstein-1d-w-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "w_{1D}", \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "v_{2D}", \
        "../../data/holstein/variational/model/holstein-2d-w-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "w_{2D}", \
        "../../data/holstein/variational/model/holstein-3d-v-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "v_{3D}", \
        "../../data/holstein/variational/model/holstein-3d-w-alpha-0to12-omega-0to2-beta-inf.dat" u 1:11 w l t "w_{3D}"

load "../gnuplot-render.gpt"