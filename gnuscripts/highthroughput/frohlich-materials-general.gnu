reset session

prefix = "../../plots/highthroughput/dielectric"

set xlabel "Alpha"
set ylabel "Mobility Difference %"

# Set grid and style
set grid
set style data points
set logscale y
set xrange [0.1:12]

# Function to calculate the mobility difference
mobility_diff(file) = sprintf("< awk 'NR > 1 {if ($5 == \"single\") {single[$2] = $21; alpha[$2] = $11} else if ($5 == \"multi\") {multi[$2] = $21}} END {for (material in single) {if (material in multi) {print alpha[material], -100*(single[material] - multi[material])/single[material]}}}' %s", file)

# Plot command
plot mobility_diff('../../data/frohlich/variational/general/frohlich-materials-general-conduction.txt') u 1:2 w points pt 7 ps 1 notitle

load "../gnuplot-render.gpt"