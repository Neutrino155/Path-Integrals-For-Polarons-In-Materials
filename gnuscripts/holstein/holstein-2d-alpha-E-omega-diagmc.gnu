reset session

prefix = "../../plots/holstein/holstein-2d-alpha-energy-adiabaticity-diagmc"
set linetype 1 pi -1 pt 7 lc rgb '#377EB8' dt solid # blue
set linetype 2 pi -1 pt 1 lc rgb '#E41A1C' dt solid # red
set linetype 3 pi -1 pt 2 lc rgb '#4DAF4A' dt solid # green
set linetype 4 pi -1 pt 4 lc rgb '#984EA3' dt solid # purple
set linetype 5 pi -1 pt 1 lc rgb '#FF7F00' dt solid # orange
set linetype 6 pi -1 pt 2 lc rgb '#A65628' dt solid # brown
set linetype 7 pi -1 pt 4 lc rgb '#F781BF' dt solid # pink

set style line 101 lc rgb '#000000' lt 1 lw 1
set tics nomirror out scale 0.75
set format '%g'

set style line 102 lc rgb '#808080' lt 0 lw 1
set grid back ls 102

set terminal pngcairo size 1024,768 enhanced font 'Helvetica,24'
set output prefix . '.png'
set pointsize 2.0

set key right
set multiplot

set size 1, 1
set origin 0, 0
set xlabel  "α"
set ylabel  "Energy (ħω₀)"
set yrange [-14:-4]
set ytics -14,1,-4
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_01.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_03.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:2 w l lw 2 t "γ=0.1", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:4 w l lw 2 t "γ=0.3", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:6 w l lw 2 t "γ=0.5"
        

set nokey
unset grid
unset xlabel
unset ylabel
set size 0.5, 0.5       # set size of inset
set origin 0.15, 0.17
set yrange [-6:-4]
set ytics -6,1,-4
set mytics 2
set xrange [0:5]
set xtics 0,1,5

plot    "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_01.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_03.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:2 w l lw 2 t "γ=0.1", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:4 w l lw 2 t "γ=0.3", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:6 w l lw 2 t "γ=0.5"
       
unset multiplot

set terminal pdfcairo size 3in,2in enhanced font 'Helvetica,12'
set output prefix . '-COLOUR.pdf'
set style line 101 lc rgb '#000000' lt 1 lw 1
set tics nomirror out scale 0.75
set format '%g'

set style line 102 lc rgb '#808080' lt 0 lw 1
set grid back ls 102

set pointsize 0.3
set key right

set multiplot

set size 1, 1
set origin 0, 0
set xlabel  "α"
set tics font "Helvetica, 12"
set yrange [-14:-4]
set ytics -14,1,-4
set ylabel  "Energy (ħω₀)"
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_01.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_03.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:2 w l t "γ=0.1", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:4 w l t "γ=0.3", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:6 w l t "γ=0.5"
        

set nokey
unset grid
unset xlabel
unset ylabel
set tics font "Helvetica, 10"
set size 0.5, 0.5       # set size of inset
set origin 0.15, 0.2
set yrange [-6:-4]
set ytics -6,1,-4
set mytics 2
set xrange [0:5]
set xtics 0,1,5

plot    "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_01.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_03.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:2 w l t "γ=0.1", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:4 w l t "γ=0.3", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:6 w l t "γ=0.5"
       
unset multiplot

set terminal pdfcairo size 3in,2in enhanced font 'Helvetica,12'
set output prefix . '-BW.pdf'

# Via Janert's Gnuplot in Action 2nd Edition; style3.gp
# Black and White traditional Physical-review ish.
set linetype 1 pi -1 pt 1 lc black dt solid
set linetype 2 pi -1 pt 6 lc black dt (8,6)
set linetype 3 pi -1 pt 2 lc black dt (4,3)
set linetype 4 pi -1 pt 4 lc black dt (3,6)
set linetype 5 pi -1 pt 1 lc black dt (12,5,2,5,2,5)
set linetype 6 pi -1 pt 6 lc black dt (16,8)
set linetype 7 pi -1 pt 2 lc black dt (20,6,2,6)
set linetype 8 pi -1 pt 4 lc black dt (30,10)
set pointsize 0.3
set termoption font "Helvetica,12"
set termoption fontscale 0.5
set style line 101 lc rgb '#000000' lt 1 lw 1
set tics nomirror out scale 0.75
set format '%g'

set style line 102 lc rgb '#808080' lt 0 lw 1
set grid back ls 102

set key right

set multiplot

set size 1, 1
set tics font "Helvetica, 12"
set origin 0, 0
set yrange [-14:-4]
set ytics -14,1,-4
set xlabel  "α"
set ylabel  "Energy (ħω₀)"
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_01.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_03.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:2 w l t "γ=0.1", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:4 w l t "γ=0.3", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:6 w l t "γ=0.5"
        

set nokey
unset grid
unset xlabel
unset ylabel
set tics font "Helvetica, 10"
set size 0.5, 0.5       # set size of inset
set origin 0.15, 0.2
set yrange [-6:-4]
set ytics -6,1,-4
set mytics 2
set xrange [0:5]
set xtics 0,1,5

plot    "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_01.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_03.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:2 w l t "γ=0.1", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:4 w l t "γ=0.3", \
        "../../data/holstein/variational/model/holstein-2d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:6 w l t "γ=0.5"
       
unset multiplot

set terminal dumb
set output