reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-vw-alpha-0to12-beta-inf-N-6"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "α"
set ylabel  "Variational Parameters (ω₀), N = 6"
set logscale y 2
set for [i=0:12] ytics (sprintf("2^{%d}", i) 2**i) offset 0.7,0
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:17 w l t "v₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:17 w l t "w₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:18 w l t "v₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:18 w l t "w₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:19 w l t "v₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:19 w l t "w₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:20 w l t "v₄", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:20 w l t "w₄", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:21 w l t "v₅", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:21 w l t "w₅", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:22 w l t "v₆", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:22 w l t "w₆", \

load "../gnuplot-render.gpt"
