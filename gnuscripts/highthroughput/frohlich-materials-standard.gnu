reset session

prefix = "../../plots/highthroughput/dielectric"
set key top right
set grid
set xlabel  "Total"
set ylabel  "alpha"
set logscale xy

plot    "../../data/frohlich/variational/standard/frohlich-materials-standard-conduction.txt" u 9:(-($16-$14)/$16) notitle, \

load "../gnuplot-render.gpt"