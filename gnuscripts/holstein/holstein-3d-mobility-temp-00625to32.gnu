reset session

prefix = "../../plots/holstein/holstein-3d-mobility-temp-00625to32"
set key Left right top
set grid
set xlabel  "T (ħω₀/kB)"
set ylabel  "Mobility (eω₀/m)"
set xrange [0.0625:32]
set yrange[1e-3:1e6]
set logscale y 10
set logscale x 2
unset xtics
set xtics nomirror
do for [i=-4:5] {
    set xtics add (sprintf("2^{%d}", i) 2**i)
}
unset ytics
set ytics nomirror
do for [i=-3:6] {
    set ytics add (sprintf("10^{%d}", i) 10**i)
}

plot    "../../data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-beta-003125to16.dat" u (1.0/$1):2 w l t "{/Symbol a}=0.1", \
        "../../data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-beta-003125to16.dat" u (1.0/$1):21 w l t "{/Symbol a}=2", \
        "../../data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-beta-003125to16.dat" u (1.0/$1):41 w l t "{/Symbol a}=4", \
        "../../data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-beta-003125to16.dat" u (1.0/$1):61 w l t "{/Symbol a}=6", \
        "../../data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-beta-003125to16.dat" u (1.0/$1):81 w l t "{/Symbol a}=8", \
        "../../data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-beta-003125to16.dat" u (1.0/$1):101 w l t "{/Symbol a}=10", \
        "../../data/holstein/variational/model/holstein-3d-mobility-alpha-0to12-beta-003125to16.dat" u (1.0/$1):121 w l t "{/Symbol a}=12"
     

load "../gnuplot-render.gpt"