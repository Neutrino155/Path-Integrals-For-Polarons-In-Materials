reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-vw-alpha-7-temp-00325to32-N-1"
set key outside center right
set rmargin at screen 0.8
set grid
set xlabel  "T (ħω₀/kB)"
set ylabel  "Variational Parameters (ω₀), N = 1"
set logscale yx 2
set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i) 
set for [i=0:6] ytics (sprintf("4^{%d}", i) 4**i) offset 0.7,0

plot    "../../data/frohlich/variational/model/frohlich-3d-v-alpha-0to12-temp-00625to32.dat" u 1:71 w l t "v₁", \
        "../../data/frohlich/variational/model/frohlich-3d-w-alpha-0to12-temp-00625to32.dat" u 1:71 w l t "w₁", \

load "../gnuplot-render.gpt"
