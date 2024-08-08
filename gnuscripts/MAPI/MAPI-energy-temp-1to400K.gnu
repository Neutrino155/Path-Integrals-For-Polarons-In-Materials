reset session

prefix = "../../plots/MAPI/MAPI-energy-temp-1to400K"
set key top right
set grid
set xlabel  "T (K)"
set ylabel  "Polaron Binding Energy (meV)"
set xrange [1:362]
# set xtics (1,50,100,150,200,250,300,350,400)
set logscale xy 2

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.dat" u 1:(-$4) w l t "Single mode", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.dat" u 1:(-$4) w l t "Multimode", \
       
load "../gnuplot-render.gpt"