reset session

# Set the output terminal and file
prefix = "../../plots/frohlich/frohlich-3d-imag-memory-alpha-0to12-temp-001-freq-0to30-contourf"

# Set the titles
set xlabel "Frequency (ω₀)" offset 0,0.8
set ylabel "α" offset -0.5,0
set cblabel "Imag Memory (ω₀m₀⁻¹)" offset 1.0,0

# Enable grid
set grid

set origin 0.02,-0.04
set size 0.9,1.14

set yrange [0:12]

set ytics 0,1,12 offset 0.5,0
set xrange [0:30]
set xtics 0,3,30 offset 0,0.5

set logscale cb 
set cbrange [0.001:1000]

set for [i=-3:3] cbtics (sprintf("10^{%d}", i) 10**i)

set autoscale fix

set pm3d map

# Plot the preprocessed data
splot '../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-12-beta-100-freq-0to30omega.dat' matrix nonuniform notitle

load "../gnuplot-render.gpt"
