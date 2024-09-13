reset session

# Set the output terminal and file
prefix = "../../plots/multivar/frohlich-3d-real-conductivity-alpha-0to12-temp-001-freq-0to30-N-2-contourf"

# Set the titles
set xlabel "Frequency (ω₀)" offset 0,0.8
set ylabel "α" offset -0.5,0
set cblabel "Frohlich Real Conductivity (e²ħ⁻¹)" offset 1,0

# Enable grid
set grid

set origin 0.02,-0.04
set size 0.9,1.14

set yrange [0:12]

set ytics 0,1,12 offset 0.5,0
set xrange [0:30]
set xtics 0,3,30 offset 0,0.5

set logscale cb 
set cbrange [0.00001:1]

set for [i=-5:0] cbtics (sprintf("10^{%d}", i) 10**i)

set autoscale fix

input_file = '../../data/frohlich/variational/multivar/frohlich-3d-real-conductivity-alpha-0to12-beta-100-freq-0to30-N-2.dat'
transposed_file = 'transposed_matrix_data.dat'

system(sprintf("awk '{ for (i = 1; i <= NF; i++) { a[NR, i] = $i; if (i > p) p = i } } END { for (j = 1; j <= p; j++) { for (i = 1; i <= NR; i++) printf(\"%%s \", a[i, j]); printf(\"\\n\") } }' %s > %s", input_file, transposed_file))

set pm3d map

# Plot the preprocessed data
splot transposed_file matrix nonuniform notitle

load "../gnuplot-render.gpt"

system('rm ' . transposed_file)