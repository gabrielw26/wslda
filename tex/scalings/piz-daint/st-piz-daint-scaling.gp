set term post enh eps color dashed 20
set output "st-piz-daint-scaling.eps"

set style line 1 pt 5 ps 2.0 lw 2. lt 1 lc 1   
set style line 2 pt 7 ps 2.0 lw 2. lt 1 lc 3
set style line 12 pt 2 ps 1.5 lw 4. lt 1 lc -1 dt 2 

set key top left font ",16"
set key spacing 1.5
set logscale xy
set xrange[50000:1e6]
set grid

set xlabel "BdG matrix size [=2N_xN_y]"
set ylabel "cost per diagonalizaion [node-hours]"
set title "Scaling for Piz Daint"

a=1
ideal(N) = a*N**3
coeffa = 0.52 / ideal(65536)
a=coeffa

fit [0.06e6:1.0e6] ideal(x) "st-piz-daint-scaling.txt" u 1:2 i 0 via a

ilabel = sprintf("ideal scaling: %.4g (2N_xN_y)^3", a)
plot "st-piz-daint-scaling.txt" u 1:2 i 0 title "double, ELPA (STAGE2, GPU)" w p ls 1, \
     ideal(x) title ilabel  w l ls 12

set out
set terminal X

system("convert -colorspace sRGB -density 200 st-piz-daint-scaling.eps st-piz-daint-scaling.png")
