reset session

# Set the output terminal and file
prefix = "../../plots/holstein/holstein-2d-real-conductivity-alpha-5-temp-00625to32-freq-0to30-contourf"

# Set the titles
set xlabel "Frequency (ω₀)" offset 0,0.8
set ylabel "Temperature (ħω₀/kB)" offset -0.5,0
set cblabel "Holstein Real Conductivity (e²ħ⁻¹)" offset 1,0

# Enable grid
set grid

set xtics 0,2,20 offset 0,0.5

set logscale y 2
set for [i=-4:5] ytics (sprintf("2^{%d}", i) 2**i)
set ytics offset 1,0

set origin 0.02,-0.04
set size 0.9,1.14

#set yrange [0.0625:4]
set cbrange [0:0.6]
set xrange [0:20]

set autoscale fix

set pm3d map

# Plot the preprocessed data
splot '../../data/holstein/variational/model/holstein-2d-real-conductivity-alpha-5-temp-00625to32-freq-0to30.dat' matrix nonuniform notitle

load "../gnuplot-render.gpt"
