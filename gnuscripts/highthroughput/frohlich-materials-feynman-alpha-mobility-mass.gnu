reset session

prefix = "../../plots/highthroughput/feynman/feynman-alpha-mobility-mass"
set key Left left bottom
set grid
set xlabel "α" offset 0,0.8
set ylabel "Polaron Mobility (cm²V⁻¹s⁻¹)" offset -0.5,0
set cblabel "Band Mass (mₑ)" offset 0.2,0
set origin 0.03,-0.04
set size 0.9,1.14

set logscale xy
set logscale cb
set yrange [1e-6:1e6]
set xrange[1e-3:1e3]
set cbrange [0:100]

set for [i=-6:6:2] ytics (sprintf("10^{%d}", i) 10**i) offset 1,0
set for [i=-3:3] xtics (sprintf("10^{%d}", i) 10**i) offset 0,0.5
set for [i=-1:2] cbtics (sprintf("10^{%d}", i) 10**i) offset -0.5,0

set arrow 1 at 6, graph 0 to 6, graph 1 nohead lc "red" dt 3 lw 2
set label 1 at 6, graph 1 "α=6" font ",12" tc 'red' offset 0.5,-0.7

conduction_file = "../../data/frohlich/variational/feynman/frohlich-materials-feynman-conduction.txt"
valence_file = "../../data/frohlich/variational/feynman/frohlich-materials-feynman-valence.txt"

set hidden3d
set view map 

splot   conduction_file using 11:20:9 w points pt 6 lc palette notitle, \
        valence_file using 11:19:9 w points pt 4 lc palette notitle, \


load "../gnuplot-render.gpt"

