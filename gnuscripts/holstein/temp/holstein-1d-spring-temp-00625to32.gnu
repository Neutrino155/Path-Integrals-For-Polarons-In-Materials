reset session

prefix = "../../plots/holstein/holstein-1d-spring-temp-00625to32"
set key left left bottom
set grid
set xlabel "Temperature (ħω₀/kB)" offset 0,0.5
set ylabel "Holstein Spring (m₀ω₀²)" offset 0,0
set xrange [0.0625:16]
set yrange [0.1:200]
set logscale x 2
set logscale y

set for [i=-2:2] ytics (sprintf("10^{%d}", i) 10**i)
set ytics offset 0.5,0

set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i)
set xtics offset 0,0.3

plot    "../../data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat" u 1:2 w l t "{/Symbol a}=0.1", \
        "../../data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat" u 1:21 w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat" u 1:41 w l t "{/Symbol a}=4", \
        "../../data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat" u 1:61 w l t "{/Symbol a}=6", \
        "../../data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat" u 1:81 w l t "{/Symbol a}=8", \
        "../../data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat" u 1:101 w l t "{/Symbol a}=10", \
        "../../data/holstein/variational/model/holstein-1d-spring-alpha-0to12-temp-00625to32.dat" u 1:121 w l t "{/Symbol a}=12"      

load "../gnuplot-render.gpt"