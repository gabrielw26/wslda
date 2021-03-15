set term post enh eps color dashed 20
set output "st-summit-scaling.eps"

set style line 1 pt 5 ps 2.0 lw 2. lt 1 lc 1       
set style line 2 pt 2 ps 1.5 lw 4. lt 1 lc 3  

set key top left font ",16"
set key spacing 1.5
set logscale xy
set xrange[50000:2e6]
set grid

set xlabel "BdG matrix size [=2N_xN_yN_z]"
set ylabel "cost per diagonalizaion [node-hours]"
set title "Scaling for Summit"

a=1
ideal(N) = a*N**3
coeffa = 0.52 / ideal(65536)
a=coeffa

fit [0.06e6:1.0e6] ideal(x) "st-summit-scaling.txt" u 1:2 via a

ilabel = sprintf("ideal scaling: %.4g (2N_xN_yN_z)^3", a)
plot "st-summit-scaling.txt" u 1:2 title "double complex, ELPA (STAGE1, GPU)" w p ls 1, \
     ideal(x) title ilabel  w l ls 2

set out
set terminal X

system("convert -colorspace sRGB -density 200 st-summit-scaling.eps st-summit-scaling.png")
