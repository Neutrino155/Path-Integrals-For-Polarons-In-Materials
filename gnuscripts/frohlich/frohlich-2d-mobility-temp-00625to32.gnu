reset session

prefix = "../../plots/frohlich/frohlich-2d-mobility-temp-00625to32"
set key Left right top
set grid
set xlabel  "Temperature (ħω₀/kB)"
set ylabel  "Frohlich Mobility (eω₀/m)"
set xrange [0.0625:32]
set yrange[5e-5:1e6]
set logscale y 10
set logscale x 2
set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i)
set for [i=-4:6:2] ytics (sprintf("10^{%d}", i) 10**i)

plot    "../../data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat" u 1:2 w l t "{/Symbol a}=0.1", \
        "../../data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat" u 1:21 w l t "{/Symbol a}=2", \
        "../../data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat" u 1:41 w l t "{/Symbol a}=4", \
        "../../data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat" u 1:61 w l t "{/Symbol a}=6", \
        "../../data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat" u 1:81 w l t "{/Symbol a}=8", \
        "../../data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat" u 1:101 w l t "{/Symbol a}=10", \
        "../../data/frohlich/variational/model/frohlich-2d-mobility-alpha-0to12-temp-00625to32.dat" u 1:121 w l t "{/Symbol a}=12", \
        # "../../data/frohlich/diagmc/frohlich_mobility_temp_0125to8_alpha_25.txt" u 2:3:($3-$4):($3+$5) with errorbars notitle, \
        # "../../data/frohlich/diagmc/frohlich_mobility_temp_0125to8_alpha_4.txt" u 2:3:($3-$4):($3+$5) with errorbars notitle, \
        # "../../data/frohlich/diagmc/frohlich_mobility_temp_0125to8_alpha_6.txt" u 2:3:($3-$4):($3+$5) with errorbars notitle

load "../gnuplot-render.gpt"