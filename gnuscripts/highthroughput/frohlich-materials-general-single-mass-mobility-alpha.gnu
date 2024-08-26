reset session

prefix = "../../plots/highthroughput/general/general-single-mass-mobility-alpha"
set key Left left bottom
set grid
set xlabel "Band Mass (mₑ)" offset 0,0.8
set cblabel "α" offset 0.2,0
set ylabel "Polaron Mobility (cm²V⁻¹s⁻¹)" offset -0.5,0
set origin 0.03,-0.04
set size 0.9,1.14

set logscale xy
set logscale cb
set yrange [1e-6:1e6]
set cbrange[1e-2:1e2]
set xrange [1e-3:1e3]

set for [i=-6:6:2] ytics (sprintf("10^{%d}", i) 10**i) offset 1,0 
set for [i=-2:2] cbtics (sprintf("10^{%d}", i) 10**i) offset -0.5,0
set for [i=-3:3] xtics (sprintf("10^{%d}", i) 10**i) offset 0,0.5

conduction_file = "../../data/frohlich/variational/general/frohlich-materials-general-single-conduction.txt"
valence_file = "../../data/frohlich/variational/general/frohlich-materials-general-single-valence.txt"

set hidden3d
set view map 

splot   conduction_file using 9:21:11 w points pt 6 lc palette notitle, \
        valence_file using 9:21:11 w points pt 4 lc palette notitle, \


load "../gnuplot-render.gpt"

