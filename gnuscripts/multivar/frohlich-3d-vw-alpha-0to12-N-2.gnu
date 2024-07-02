reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-vw-alpha-0to12-beta-inf-N-2"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "α"
set ylabel  "Variational Parameters (ω₀), N = 2"
set logscale y 2
set for [i=0:6] ytics (sprintf("2^{%d}", i) 2**i) offset 0.7,0
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:3 w l t "v₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:3 w l t "w₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:4 w l t "v₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:4 w l t "w₂", \
        "../../data/frohlich/variational/multivar/abe-okamoto-alpha-0.5to5-beta-inf-N-2.dat" u 1:3 pt 7 dt 1 t "v₁", \
        "../../data/frohlich/variational/multivar/abe-okamoto-alpha-0.5to5-beta-inf-N-2.dat" u 1:2 pt 7 dt 1 t "w₁", \
        "../../data/frohlich/variational/multivar/abe-okamoto-alpha-0.5to5-beta-inf-N-2.dat" u 1:5 pt 7 dt 1 t "v₂", \
        "../../data/frohlich/variational/multivar/abe-okamoto-alpha-0.5to5-beta-inf-N-2.dat" u 1:4 pt 7 dt 1 t "w₂", \

load "../gnuplot-render.gpt"
