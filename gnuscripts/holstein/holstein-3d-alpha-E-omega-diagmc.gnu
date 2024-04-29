reset session

prefix = "../../plots/holstein/holstein-3d-alpha-energy-adiabaticity-diagmc"
set key right
set grid
set xlabel  "α"
set ylabel  "Energy (J)"
set xrange [0:12]
set xtics 0,1,12

plot    "../../data/holstein/variational/model/holstein-3d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:2 w l t "γ=0.1", \
        "../../data/holstein/variational/model/holstein-3d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:4 w l t "γ=0.3", \
        "../../data/holstein/variational/model/holstein-3d-E-alpha-0to12-omega-0to2-beta-inf.dat" u 1:6 w l t "γ=0.5", \
        "../../data/holstein/diagmc/holstein_3d_parabolic_energy_alpha_0to5_gamma_01.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_3d_parabolic_energy_alpha_0to5_gamma_03.txt" u 1:2:3 with yerrorbars notitle, \
        "../../data/holstein/diagmc/holstein_3d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars notitle
       
load "../gnuplot-render.gpt"