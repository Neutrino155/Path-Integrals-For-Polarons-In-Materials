reset session

prefix = "../../plots/MAPI/MAPI-imag-memory-temp-1to400K"
set key right
set grid
set xlabel  "Temperature (K)"
set ylabel  "Imaginary Memory Imχ (m₀/ω₀)"
set xrange [1:362]
# set xtics (1,50,100,150,200,250,300,350,400)
set logscale xy 2

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.dat" u 1:6 w l t "Single mode", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.dat" u 1:6 w l t "Multimode"
       
load "../gnuplot-render.gpt"