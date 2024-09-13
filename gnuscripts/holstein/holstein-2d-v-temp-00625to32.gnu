reset session

prefix = "../../plots/holstein/holstein-2d-v-temp-00625to32"
set key right bottom
set grid
set xlabel  "T (ħω₀/kB)"
set ylabel  "Holstein v (ω₀)"
set xrange [0.0625:32]
set yrange [1:40]
set logscale xy 2

set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i)

plot    "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat" u 1:2 w l t "{/Symbol a}=0.1", \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat" u 1:21 w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat" u 1:41 w l t "{/Symbol a}=4", \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat" u 1:61 w l t "{/Symbol a}=6", \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat" u 1:81 w l t "{/Symbol a}=8", \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat" u 1:101 w l t "{/Symbol a}=10", \
        "../../data/holstein/variational/model/holstein-2d-v-alpha-0to12-temp-00625to32.dat" u 1:121 w l t "{/Symbol a}=12"      

load "../gnuplot-render.gpt"
