reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-vw-alpha-7-temp-00325to32-N-4"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "T (ħω₀/kB)"
set ylabel  "Variational Parameters (ω₀), N = 4"
set logscale yx 2
set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i) 
set for [i=0:7] ytics (sprintf("4^{%d}", i) 4**i) offset 0.7,0

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-7-temp-00625to32-N-4.dat" u 1:2 w l t "v₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-7-temp-00625to32-N-4.dat" u 1:2 w l t "w₁", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-7-temp-00625to32-N-4.dat" u 1:3 w l t "v₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-7-temp-00625to32-N-4.dat" u 1:3 w l t "w₂", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-7-temp-00625to32-N-4.dat" u 1:4 w l t "v₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-7-temp-00625to32-N-4.dat" u 1:4 w l t "w₃", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-v-alpha-7-temp-00625to32-N-4.dat" u 1:5 w l t "v₄", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-w-alpha-7-temp-00625to32-N-4.dat" u 1:5 w l t "w₄", \

load "../gnuplot-render.gpt"
