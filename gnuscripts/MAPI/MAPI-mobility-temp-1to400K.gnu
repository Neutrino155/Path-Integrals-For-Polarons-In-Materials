reset session

prefix = "../../plots/MAPI/MAPI-mobility-temp-1to400K"
set key right
set grid
set xlabel  "T (K)"
set ylabel  "Mobility μ (cm²V⁻¹s⁻¹)"
set xrange [1:400]
set xtics (1,50,100,150,200,250,300,350,400)
set logscale y 2
set yrange [64:1028]

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.dat" u 1:7 w l t "Single mode", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.dat" u 1:7 w l t "Multimode"
       
load "../gnuplot-render.gpt"