reset session

# Set the output terminal and file
prefix = "../../plots/holstein/holstein-3d-w-temp-00625to32-contourf"

# Set the titles
set xlabel "α" offset 0,0.8
set ylabel "Temperature (ħω₀/kB)" offset -0.5,0
set cblabel "Holstein w (ω₀)" offset 1.5,0

# Enable grid
set grid

set xrange [0:12]
set xtics 0,1,12 offset 0,0.5
set logscale y 2
set logscale cb 2
set cbrange [1:256]

set origin 0.02,-0.04
set size 0.9,1.14

set for [i=-4:5] ytics (sprintf("2^{%d}", i) 2**i)
set ytics offset 1,0

set for [i=0:8] cbtics (sprintf("2^{%d}", i) 2**i)

set autoscale fix

set pm3d map

# Plot the preprocessed data
splot '../../data/holstein/variational/model/holstein-3d-w-alpha-0to12-temp-00625to32.dat' matrix nonuniform notitle

load "../gnuplot-render.gpt"
