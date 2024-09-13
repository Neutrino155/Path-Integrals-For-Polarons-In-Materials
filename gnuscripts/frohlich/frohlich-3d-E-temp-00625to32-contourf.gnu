reset session

# Set the output terminal and file
prefix = "../../plots/frohlich/frohlich-3d-energy-temp-00625to32-contourf"

# Set the titles
set xlabel "α" offset 0,0.8
set ylabel "Temperature (ħω₀/kB)" offset -0.5,0
set cblabel "Frohlich Free Energy |F| (ħω₀)" offset 0.5,0

# Enable grid
set grid

set origin 0.02,-0.04
set size 0.9,1.14

set xrange [0:12]
set yrange [0:32]
set cbrange [2**-3:2**8]
set xtics 0,1,12 offset 0,0.5
set logscale y 2
set logscale cb 2

unset ytics
set ytics nomirror
set for [i=-4:5] ytics (sprintf("2^{%d}", i) 2**i)
set ytics offset 1,0

unset cbtics
set cbtics nomirror
set for [i=-3:8] cbtics (sprintf("2^{%d}", i) 2**i)
set cbtics offset -0.5,0

set autoscale fix

set pm3d map

# Plot the preprocessed data
splot '../../data/frohlich/variational/model/frohlich-3d-E-alpha-0to12-temp-00625to32.dat' matrix nonuniform notitle

load "../gnuplot-render.gpt"
