reset session

prefix = "../../plots/frohlich/frohlich-3d-mobility-alpha-6-temp-00625to32"
set key Left right top
set grid
set xlabel  "Temperature (ħω₀/kB)"
set ylabel  "Frohlich Mobility (eω₀/m)"
set xrange [0.0625:32]
set yrange [10**-3:10**4]
set logscale y 10
set logscale x 2
set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i)
set for [i=-3:4] ytics (sprintf("10^{%d}", i) 10**i)

set label 1 "α = 6.0" at 2**-2,3*10**3

plot    "../../data/frohlich/variational/model/frohlich-3d-mobility-alpha-0to12-temp-00625to32.dat" u 1:(1/$1) w l dt 2 t "MIR", \
        "../../data/frohlich/variational/model/frohlich-3d-mobility-alpha-0to12-temp-00625to32.dat" u 1:61 w l dt 1 t "variational", \
        "../../data/frohlich/diagmc/frohlich_mobility_temp_0125to8_alpha_6.txt" u 2:3:($3-$4):($3+$5) with errorbars pt 7 dt 1 t "diagmc", \

load "../gnuplot-render.gpt"