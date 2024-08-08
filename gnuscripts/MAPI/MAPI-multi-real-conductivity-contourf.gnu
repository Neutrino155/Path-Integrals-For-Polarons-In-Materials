reset session

# Set the output terminal and file
prefix = "../../plots/MAPI/MAPI-multi-real-conductivity-contourf"

# Set the titles
set xlabel "Frequency (THz2π)" offset 0,0.8
set ylabel "Temperature (K)" offset -0.5,0
set cblabel "Real Conductivity (μS)" offset 1,0

# Enable grid
set grid

set xtics 0,4,48 offset 0,0.5
set ytics 0,50,400 offset 0.5,0

set origin 0.02,-0.04
set size 0.9,1.14

set cbrange [0:100]

set autoscale fix

input_file = '../../data/frohlich/variational/MAPI/frohlich-MAPI-multi-real-conductivity-temp-0to400K-freq-0to30omega.dat'
transposed_file = 'transposed_matrix_data.dat'

system(sprintf("awk '{ for (i = 1; i <= NF; i++) { a[NR, i] = $i; if (i > p) p = i } } END { for (j = 1; j <= p; j++) { for (i = 1; i <= NR; i++) printf(\"%%s \", a[i, j]); printf(\"\\n\") } }' %s > %s", input_file, transposed_file))

set pm3d map

# Plot the preprocessed data
splot transposed_file matrix nonuniform notitle

load "../gnuplot-render.gpt"

system('rm ' . transposed_file)