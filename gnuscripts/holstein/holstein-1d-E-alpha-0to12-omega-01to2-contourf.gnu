reset session

prefix = "../../plots/holstein/holstein-1d-energy-alpha-0to12-omega-01to2-contourf"

# Set the titles
set xlabel "α" offset 0,0.8
set ylabel "Phonon Frequency (ω₀)" offset -0.5,0
set cblabel "Holstein Energy (ħω₀)" offset -0.2,0

# Enable grid
set grid

set origin 0.02,-0.04
set size 0.9,1.14

set xrange [0:12]
set yrange [0.1:2]
set cbrange [*:*]
set xtics 0,1,12 offset 0,0.5
set ytics 0.0,0.2,2 offset 1,0

input_file = '../../data/holstein/variational/model/holstein-1d-energy-alpha-0to12-beta-inf.dat'
transposed_file = 'transposed_matrix_data.dat'

system(sprintf("awk '{ for (i = 1; i <= NF; i++) { a[NR, i] = $i; if (i > p) p = i } } END { for (j = 1; j <= p; j++) { for (i = 1; i <= NR; i++) printf(\"%%s \", a[i, j]); printf(\"\\n\") } }' %s > %s", input_file, transposed_file))

set pm3d map

# Plot the preprocessed data
splot transposed_file matrix nonuniform notitle

load "../gnuplot-render.gpt"

system('rm ' . transposed_file)