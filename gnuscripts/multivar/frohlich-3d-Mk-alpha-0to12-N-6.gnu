reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-Mk-alpha-0to12-beta-inf-N-6"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "α"
set ylabel  "Mass (m₀) & Spring (m₀ω₀²), N = 6"
set logscale y 4
set for [i=-3:12] ytics (sprintf("4^{%d}", i) 4**i) offset 0.7,0
set xrange [0:12]
set yrange [4**-3:*]
set xtics 0,1,12

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:17 w l t "M₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:17 w l t "κ₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:18 w l t "M₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:18 w l t "κ₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:19 w l t "M₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:19 w l t "κ₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:20 w l t "M₄", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:20 w l t "κ₄", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:21 w l t "M₅", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:21 w l t "κ₅", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:22 w l t "M₆", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:22 w l t "κ₆", \

load "../gnuplot-render.gpt"
