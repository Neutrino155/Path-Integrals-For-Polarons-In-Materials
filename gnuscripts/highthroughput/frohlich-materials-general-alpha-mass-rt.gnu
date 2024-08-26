reset session

prefix = "../../plots/highthroughput/general/general-alpha-mass-rt"
set key Left left top
set grid
set xlabel  "α"
set ylabel  "300K Polaron Mass (mₑ)"

set logscale xy
set yrange [1e-2:1e10]
set xrange [1e-3:1e3]

set for [i=-2:10:2] ytics (sprintf("10^{%d}", i) 10**i)
set for [i=-3:3] xtics (sprintf("10^{%d}", i) 10**i)

set arrow 1 at 6, graph 0 to 6, graph 1 nohead lc "red" dt 3 lw 2
set label 1 at 6, graph 1 "α=6" font ",12" tc 'red' offset 0.5,-0.7

conduction_file = "../../data/frohlich/variational/general/frohlich-materials-general-single-conduction.txt"
valence_file = "../../data/frohlich/variational/general/frohlich-materials-general-single-valence.txt"

system(sprintf("awk '$4 ~ /%s/ {print}' %s > %s", '15', conduction_file, 'conduction_group_15_output.txt'))
system(sprintf("awk '$4 ~ /%s/ {print}' %s > %s", '16', conduction_file, 'conduction_group_16_output.txt'))
system(sprintf("awk '$4 ~ /%s/ {print}' %s > %s", '17', conduction_file, 'conduction_group_17_output.txt'))
system(sprintf("awk '$4 !~ /%s/ {print}' %s > %s", '15', conduction_file, 'conduction_no_group_15_output.txt'))
system(sprintf("awk '$4 !~ /%s/ {print}' %s > %s", '16', 'conduction_no_group_15_output.txt', 'conduction_no_group_16_output.txt'))
system(sprintf("awk '$4 !~ /%s/ {print}' %s > %s", '17', 'conduction_no_group_16_output.txt', 'conduction_other_groups_output.txt'))

system(sprintf("awk '$4 ~ /%s/ {print}' %s > %s", '15', valence_file, 'valence_group_15_output.txt'))
system(sprintf("awk '$4 ~ /%s/ {print}' %s > %s", '16', valence_file, 'valence_group_16_output.txt'))
system(sprintf("awk '$4 ~ /%s/ {print}' %s > %s", '17', valence_file, 'valence_group_17_output.txt'))
system(sprintf("awk '$4 !~ /%s/ {print}' %s > %s", '15', valence_file, 'valence_no_group_15_output.txt'))
system(sprintf("awk '$4 !~ /%s/ {print}' %s > %s", '16', 'valence_no_group_15_output.txt', 'valence_no_group_16_output.txt'))
system(sprintf("awk '$4 !~ /%s/ {print}' %s > %s", '17', 'valence_no_group_16_output.txt', 'valence_other_groups_output.txt'))


plot    "conduction_other_groups_output.txt" using 11:($14**2/$15**2*$9) w points pt 6 lc rgb '#984EA3' t 'Other', \
        "conduction_group_15_output.txt" using 11:($14**2/$15**2*$9) w points pt 6 lc rgb '#4DAF4A' t '15', \
        "conduction_group_16_output.txt" using 11:($14**2/$15**2*$9) w points pt 6 lc rgb '#FF7F00' t '16', \
        "conduction_group_17_output.txt" using 11:($14**2/$15**2*$9) w points pt 6 lc rgb '#377EB8' t '17', \
        "valence_other_groups_output.txt" using 11:($14**2/$15**2*$9) w points pt 4 lc rgb '#984EA3' notitle, \
        "valence_group_15_output.txt" using 11:($14**2/$15**2*$9) w points pt 4 lc rgb '#4DAF4A' notitle, \
        "valence_group_16_output.txt" using 11:($14**2/$15**2*$9) w points pt 4 lc rgb '#FF7F00' notitle, \
        "valence_group_17_output.txt" using 11:($14**2/$15**2*$9) w points pt 4 lc rgb '#377EB8' notitle, \

load "../gnuplot-render.gpt"

system('rm *_output.txt')