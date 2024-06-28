reset session

prefix = "../../plots/frohlich/frohlich-3d-imag-memory-freq-3d"
unset key
set xlabel  "alpha"
set ylabel  "Frequency (ω₀)"
set zlabel  "Frohlich Imag Memory (m/ω₀)"
set grid

set logscale z
set zrange [.1:*]

stats '../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-1to7-beta-100-freq-0to30.dat' using 1 nooutput
rows = STATS_records

stats '../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-1to7-beta-100-freq-0to30.dat' matrix nooutput
cols = STATS_records/rows

set view 60,40
set xyplane 0

splot for [i=1:cols] '../../data/frohlich/variational/model/frohlich-3d-imag-memory-alpha-1to7-beta-100-freq-0to30.dat' matrix every ::i::i lt 1 lw 2 with lines notitle

load "../gnuplot-render.gpt"