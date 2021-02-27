set term post enh eps color dashed 20
set output "td-summit-scaling-1.eps"

set style line 1 pt 5 ps 2.0 lw 2. lt 1 lc 1       
set style line 2 pt 2 ps 1.5 lw 4. lt 1 lc 3  

set key top left font ",16"
set key spacing 1.5
set logscale xy
set xrange[10000:1e6]
set grid

set xlabel "Number of lattice points N [=N_xN_yN_z]"
set ylabel "cost per trajectory of lenght te_F=1 \n[node-hours]" offset 3,0
set title "Summit, td-wslda-3d, SLDA, spin symmetric system"

a=1
ideal(N) = a*N**2*log(N)
coeffa = 0.0346   / ideal(32*32*32)
a=coeffa

fit [0.06e6:0.8e6] ideal(x) "td-summit-scaling-1.txt" u ($1*$2*$3):4 via a

ilabel = sprintf("ideal scaling: %.4g N^2log(N)", a)
plot "td-summit-scaling-1.txt" u ($1*$2*$3):4 title "measurment" w p ls 1, \
     ideal(x) title ilabel  w l ls 2

set out
set terminal X

system("convert -colorspace sRGB -density 200 td-summit-scaling-1.eps td-summit-scaling-1.png")
