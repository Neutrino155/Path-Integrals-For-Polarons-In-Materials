reset session

prefix = "../../plots/multivar/frohlich-3d-multivariate-energy-alpha-7-temp-00078125to32-N-1to5"
set key top right
set grid
set xlabel  "T (ħω₀/kB)"
set ylabel  "Absolute Difference with N = 1"
set logscale x 4
set for [i=-4:3] xtics (sprintf("4^{%d}", i) 4**i) 

plot    "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-7-temp-00078125to32-N-1to5.dat" u 1:($2-$6) w l t "N=5", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-7-temp-00078125to32-N-1to5.dat" u 1:($2-$5) w l t "N=4", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-7-temp-00078125to32-N-1to5.dat" u 1:($2-$4) w l t "N=3", \
        "../../data/frohlich/variational/multivar/frohlich-3d-multivariate-energy-alpha-7-temp-00078125to32-N-1to5.dat" u 1:($2-$3) w l t "N=2", \

load "../gnuplot-render.gpt"
