reset session

prefix = "../../plots/MAPI/MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi"

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

set grid
set xlabel  "Frequency (THz2π)"
set ylabel  "Real Conductivity σ (μS)"
set xrange [0:100]
set mxtics 2
set yrange [-1:160]
set size 1, 1
set origin 0, 0

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:2 w l t "0.48 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:4 w l t "100 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:6 w l t "200 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:8 w l t "300 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:10 w l t "400 K"
       
set nokey
unset grid
unset xlabel
unset ylabel
set size 0.45, 0.5       # set size of inset
set origin 0.25, 0.45
set xrange [0:5]
set yrange [-1:160]

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:2 w l t "0.48 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:4 w l t "100 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:6 w l t "200 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:8 w l t "300 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:10 w l t "400 K"

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

set grid
set xlabel  "Frequency (THz2π)"
set ylabel  "Real Conductivity σ (μS)"
set xrange [0:100]
set mxtics 2
set yrange [-1:160]
set size 1, 1
set origin 0, 0

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:2 w l t "0.48 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:4 w l t "100 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:6 w l t "200 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:8 w l t "300 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:10 w l t "400 K"
       
set nokey
unset grid
unset xlabel
unset ylabel
set size 0.45, 0.5       # set size of inset
set origin 0.25, 0.45
set tics font "Helvetica, 10"
set xrange [0:5]
set yrange [-1:160]

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:2 w l t "0.48 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:4 w l t "100 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:6 w l t "200 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:8 w l t "300 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:10 w l t "400 K"

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

set grid
set xlabel  "Frequency (THz2π)"
set ylabel  "Real Conductivity σ (μS)"
set xrange [0:100]
set mxtics 2
set yrange [-1:160]
set size 1, 1
set origin 0, 0

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:2 w l t "0.48 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:4 w l t "100 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:6 w l t "200 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:8 w l t "300 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to100THz2pi.dat" u 1:10 w l t "400 K"
       
set nokey
unset grid
unset xlabel
unset ylabel
set size 0.45, 0.5       # set size of inset
set origin 0.25, 0.45
set tics font "Helvetica, 10"
set xrange [0:5]
set yrange [-1:160]

plot    "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:2 w l t "0.48 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:4 w l t "100 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:6 w l t "200 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:8 w l t "300 K", \
        "../../data/frohlich/variational/MAPI/frohlich-MAPI-single-real-conductivity-temp-0.48to400K-freq-0to10omega.dat" u 1:10 w l t "400 K"

unset multiplot

set terminal dumb
set output