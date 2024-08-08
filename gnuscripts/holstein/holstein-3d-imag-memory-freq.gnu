reset session

prefix = "../../plots/holstein/holstein-3d-imag-memory-freq"
set key Left right top
set xlabel  "Frequency (ω₀)"
set ylabel  "Holstein Imag Memory (ω₀/m₀)"
set yrange [0:10]

input_file = '../../data/holstein/variational/model/holstein-3d-imag-memory-alpha-0to12-beta-100-freq-0to30.dat'
transposed_file = 'transposed_matrix_data.dat'

system(sprintf("awk '{ for (i = 1; i <= NF; i++) { a[NR, i] = $i; if (i > p) p = i } } END { for (j = 1; j <= p; j++) { for (i = 1; i <= NR; i++) printf(\"%%s \", a[i, j]); printf(\"\\n\") } }' %s > %s", input_file, transposed_file))

plot    transposed_file u 1:21 w l t "{/Symbol a}=1", \
        transposed_file u 1:41 w l t "{/Symbol a}=2", \
        transposed_file u 1:61 w l t "{/Symbol a}=3", \
        transposed_file u 1:81 w l t "{/Symbol a}=4", \
        transposed_file u 1:101 w l t "{/Symbol a}=5", \
        transposed_file u 1:121 w l t "{/Symbol a}=6", \
        transposed_file u 1:141 w l t "{/Symbol a}=7"

load "../gnuplot-render.gpt"

system('rm ' . transposed_file)