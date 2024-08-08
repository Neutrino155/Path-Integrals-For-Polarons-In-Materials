reset session

prefix = "../../plots/holstein/holstein-1d-energy-temp-00625to32"
set key right bottom opaque
set grid
set xlabel "Temperature (ħω₀/kB)" offset 0,0.5
set ylabel "Negative Holstein Energy (-ħω₀)" offset 0.5,0
set xrange [0.0625:16]
set logscale xy 2

set for [i=-1:5] ytics (sprintf("2^{%d}", i) 2**i)
set ytics offset 0.5,0

set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i)
set xtics offset 0,0.3

plot    "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat" u 1:(-$2) w l t "{/Symbol a}=0.1", \
        "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat" u 1:(-$21) w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat" u 1:(-$41) w l t "{/Symbol a}=4", \
        "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat" u 1:(-$61) w l t "{/Symbol a}=6", \
        "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat" u 1:(-$81) w l t "{/Symbol a}=8", \
        "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat" u 1:(-$101) w l t "{/Symbol a}=10", \
        "../../data/holstein/variational/model/holstein-1d-E-alpha-0to12-temp-00625to32.dat" u 1:(-$121) w l t "{/Symbol a}=12"      

load "../gnuplot-render.gpt"