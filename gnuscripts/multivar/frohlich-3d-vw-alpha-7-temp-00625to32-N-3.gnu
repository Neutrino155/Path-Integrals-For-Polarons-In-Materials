reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-vw-alpha-7-temp-00325to32-N-3"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "T (ħω₀/kB)"
set ylabel  "Variational Parameters (ω₀), N = 3"
set logscale yx 2
set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i) 
set for [i=0:6] ytics (sprintf("4^{%d}", i) 4**i) offset 0.7,0

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v1-alpha-0to12-temp-00625to32-N-3.dat" u 1:71 w l t "v₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w1-alpha-0to12-temp-00625to32-N-3.dat" u 1:71 w l t "w₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v2-alpha-0to12-temp-00625to32-N-3.dat" u 1:71 w l t "v₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w2-alpha-0to12-temp-00625to32-N-3.dat" u 1:71 w l t "w₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v3-alpha-0to12-temp-00625to32-N-3.dat" u 1:71 w l t "v₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w3-alpha-0to12-temp-00625to32-N-3.dat" u 1:71 w l t "w₃", \

load "../gnuplot-render.gpt"
