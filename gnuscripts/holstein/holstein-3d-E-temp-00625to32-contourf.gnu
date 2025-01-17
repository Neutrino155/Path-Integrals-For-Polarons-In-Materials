reset session

# Set the output terminal and file
prefix = "../../plots/holstein/holstein-3d-energy-temp-00625to32-contourf"

# Set the titles
set xlabel "α" offset 0,0.8
set ylabel "Temperature (ħω₀/kB)" offset -0.5,0
set cblabel "Holstein Energy (ħω₀)" offset -0.2,0

# Enable grid
set grid

set origin 0.02,-0.04
set size 0.9,1.14

set xrange [0:12]
set yrange [0:32]
set cbrange [-120:0]
set xtics 0,1,12 offset 0,0.5
set ytics 0,4,32 offset 1,0
#set logscale y 2

#unset ytics
#set ytics nomirror
#set for [i=-4:5] ytics (sprintf("2^{%d}", i) 2**i)
#set ytics offset 1,0

set autoscale fix

set pm3d map

# Plot the preprocessed data
splot '../../data/holstein/variational/model/holstein-3d-E-alpha-0to12-temp-00625to32.dat' matrix nonuniform notitle

load "../gnuplot-render.gpt"
