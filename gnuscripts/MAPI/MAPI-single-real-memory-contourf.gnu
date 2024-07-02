reset session

# Set the output terminal and file
prefix = "../../plots/MAPI/MAPI-single-real-memory-contourf"

# Set the titles
set xlabel "Frequency (THz2π)" offset 0,0.8
set ylabel "Temperature (K)" offset -0.5,0
set cblabel "Real Memory (ω₀/m₀)" offset 1,0

# Enable grid
set grid

set xtics 0,2,22 offset 0,0.5
set ytics 0,50,400 offset 0.5,0

set origin 0.02,-0.04
set size 0.9,1.14

set cbrange [-7.2:3]

set autoscale fix

set pm3d map

# Plot the preprocessed data
splot '../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-memory-temp-0to400K-freq-0to22.5THz2pi.dat' matrix nonuniform notitle

load "../gnuplot-render.gpt"
