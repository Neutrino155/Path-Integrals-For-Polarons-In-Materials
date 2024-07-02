reset session

prefix = "../../plots/MAPI/MAPI-multi-imag-conductivity-temp-0.48to400K-freq-0to10omega"
set key top right
set grid
set xlabel  "Frequency (THz2π)"
set ylabel  "Imag Conductivity σ (μS)"
set yrange [-40:180]
set xrange [0:22.5]
set mxtics 2

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:(-$2) w l t "0.48 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:(-$4) w l t "100 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:(-$6) w l t "200 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:(-$8) w l t "300 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-imag-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:(-$10) w l t "400 K"
       
load "../gnuplot-render.gpt"