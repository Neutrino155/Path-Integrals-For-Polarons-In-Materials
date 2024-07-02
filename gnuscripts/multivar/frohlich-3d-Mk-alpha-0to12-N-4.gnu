reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-Mk-alpha-0to12-beta-inf-N-4"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "α"
set ylabel  "Mass (m₀) & Spring (m₀ω₀²), N = 4"
set logscale y 4
set for [i=-3:10] ytics (sprintf("4^{%d}", i) 4**i) offset 0.7,0
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:8 w l t "M₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:8 w l t "κ₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:9 w l t "M₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:9 w l t "κ₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:10 w l t "M₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:10 w l t "κ₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-mass-alpha-0to12-beta-inf-N-1to6.dat" u 1:11 w l t "M₄", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-spring-alpha-0to12-beta-inf-N-1to6.dat" u 1:11 w l t "κ₄", \

load "../gnuplot-render.gpt"
