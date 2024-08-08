reset session

prefix = "../../plots/highthroughput/dielectric"
set key top right
set grid
set xlabel  "Total"
set ylabel  "alpha"

plot    "../../data/frohlich/variational/feynman/frohlich-materials-feynman-conduction.txt" u 9:(($16-$15)/$16) notitle, \

load "../gnuplot-render.gpt"