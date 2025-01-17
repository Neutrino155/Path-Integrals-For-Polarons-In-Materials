reset session

prefix = "../../plots/holstein/holstein-2d-energy-temp-00625to32"
set key Left left bottom
set grid
set xlabel  "T (ħω₀/kB)"
set ylabel  "Holstein Energy (ħω₀)"
set xrange [0:32]
set xtics 0,4,32
set yrange [*:*]
#set logscale x 2
#unset xtics
#set xtics nomirror
#do for [i=-4:5] {
#    set xtics add (sprintf("2^{%d}", i) 2**i)
#}

plot    "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat" u 1:2 w l t "{/Symbol a}=0.1", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat" u 1:21 w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat" u 1:41 w l t "{/Symbol a}=4", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat" u 1:61 w l t "{/Symbol a}=6", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat" u 1:81 w l t "{/Symbol a}=8", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat" u 1:101 w l t "{/Symbol a}=10", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-temp-00625to32.dat" u 1:121 w l t "{/Symbol a}=12"      

load "../gnuplot-render.gpt"
