reset session

prefix = "../../plots/MAPI/MAPI-mobility-temp-1to400K"
set key right
set grid
set xlabel  "Temperature (K)"
set ylabel  "Mobility μ (cm²V⁻¹s⁻¹)"
set xtics (1,50,100,150,200,250,300,350,400)
set yrange [100:400]
#set logscale y 10
#set for [i=-1:6] ytics (sprintf("10^{%d}", i) 10**i)


plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-temp-1to400K.dat" u 1:10 w l t "Single mode", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-temp-1to400K.dat" u 1:10 w l t "Multimode"
       
load "../gnuplot-render.gpt"