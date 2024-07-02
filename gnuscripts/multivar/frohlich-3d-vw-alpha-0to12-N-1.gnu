reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-vw-alpha-0to12-beta-inf-N-1"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "α"
set ylabel  "Variational Parameters (ω₀), N = 1"
set logscale y 2
set for [i=0:5] ytics (sprintf("2^{%d}", i) 2**i) offset 0.7,0
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-0to12-beta-inf-N-1to6.dat" u 1:2 w l t "v₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-0to12-beta-inf-N-1to6.dat" u 1:2 w l t "w₁", \

load "../gnuplot-render.gpt"
