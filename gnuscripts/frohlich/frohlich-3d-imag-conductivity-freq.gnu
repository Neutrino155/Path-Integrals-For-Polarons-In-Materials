reset session

prefix = "../../plots/frohlich/frohlich-3d-imag-conductivity-freq"
set key Left right top
set xlabel  "Frequency (ω₀)"
set ylabel  "Frohlich Imag Optical Absorption (ϵ₀cn)⁻¹"
set yrange [-0.5:0.9]

plot    "../../data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-1to7-beta-100-freq-0to30.dat" u 1:2 w l t "{/Symbol a}=1", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-1to7-beta-100-freq-0to30.dat" u 1:3 w l t "{/Symbol a}=2", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-1to7-beta-100-freq-0to30.dat" u 1:5 w l t "{/Symbol a}=3", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-1to7-beta-100-freq-0to30.dat" u 1:6 w l t "{/Symbol a}=4", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-1to7-beta-100-freq-0to30.dat" u 1:7 w l t "{/Symbol a}=5", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-1to7-beta-100-freq-0to30.dat" u 1:8 w l t "{/Symbol a}=6", \
        "../../data/frohlich/variational/model/frohlich-3d-imag-conductivity-alpha-1to7-beta-100-freq-0to30.dat" u 1:9 w l t "{/Symbol a}=7"

load "../gnuplot-render.gpt"