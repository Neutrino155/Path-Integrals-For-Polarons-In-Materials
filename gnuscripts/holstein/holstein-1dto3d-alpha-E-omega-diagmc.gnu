reset session

prefix = "../../plots/holstein/holstein-1dto3d-alpha-energy-omega-diagmc"
set key Left left bottom
set xlabel  "α"
set ylabel  "Energy (ħω₀)"
set xrange [0:5]
set xtics 0,1,5

plot    "../../data/holstein/diagmc/holstein_1d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars t "1D PB", \
        "../../data/holstein/diagmc/holstein_2d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars t "2D PB", \
        "../../data/holstein/diagmc/holstein_3d_parabolic_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars t "3D PB", \
        "../../data/holstein/diagmc/holstein_1d_tightbinding_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars t "1D TB", \
        "../../data/holstein/diagmc/holstein_2d_tightbinding_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars t "2D TB", \
        "../../data/holstein/diagmc/holstein_3d_tightbinding_energy_alpha_0to5_gamma_05.txt" u 1:2:3 with yerrorbars t "3D TB", \
        "../../data/holstein/variational/model/holstein-1d-energy-alpha-0to12-beta-inf.dat" u 1:6 w l lw 2 t "1D V", \
        "../../data/holstein/variational/model/holstein-2d-energy-alpha-0to12-beta-inf.dat" u 1:6 w l lw 2 t "2D V", \
        "../../data/holstein/variational/model/holstein-3d-energy-alpha-0to12-beta-inf.dat" u 1:6 w l lw 2 t "3D V"
     
load "../gnuplot-render.gpt"