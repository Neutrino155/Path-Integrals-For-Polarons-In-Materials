reset session

prefix = "../../plots/MAPI/MAPI-vw-temp-1to400K"
set key bottom right
set grid
set xlabel  "T (K)"
set ylabel  "v, w (ω₀)"
set xrange [1:362]
#set xtics (1,50,100,150,200,250,300,350,400)
set logscale xy 2


plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.dat" u 1:2 w l t "v_{Single}", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.dat" u 1:3 w l t "w_{Single}", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.dat" u 1:2 w l t "v_{Multi}", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.dat" u 1:3 w l t "w_{Multi}"
       
load "../gnuplot-render.gpt"