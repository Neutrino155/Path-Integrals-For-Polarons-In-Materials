reset session

prefix = "../../plots/frohlich/frohlich-2d-energy-temp-00625to32"
set key at 4.5,0.7 horizontal maxcols 2
set grid
set xlabel  "Temperature (ħω₀/kB)"
set ylabel  "Frohlich Free Energy |F| (ħω₀)"
set xrange [0:32]
set yrange [2**-3:2**8]
set logscale xy 2
set for [i=-4:5] xtics (sprintf("2^{%d}", i) 2**i)
set for [i=-3:8] ytics (sprintf("2^{%d}", i) 2**i)


plot    "../../data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat" u 1:($2) w l t "{/Symbol a}=0.1", \
        "../../data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat" u 1:($21) w l t "{/Symbol a}=2", \
        "../../data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat" u 1:($41) w l t "{/Symbol a}=4", \
        "../../data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat" u 1:($61) w l t "{/Symbol a}=6", \
        "../../data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat" u 1:($81) w l t "{/Symbol a}=8", \
        "../../data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat" u 1:($101) w l t "{/Symbol a}=10", \
        "../../data/frohlich/variational/model/frohlich-2d-E-alpha-0to12-temp-00625to32.dat" u 1:($121) w l t "{/Symbol a}=12"      

load "../gnuplot-render.gpt"
