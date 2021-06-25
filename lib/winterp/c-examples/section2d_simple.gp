set term png size 1024,512 font ",22"
set output "section2d_simple.png"

set style line 1 pt 7 ps 2.0 lw 2. lt 1 lc 1
set style line 2 pt 2 ps 1.5 lw 2. lt 1 lc 3 

set xrange[0:32]
plot "section2d_simple.txt" u 1:2 i 0 t "lattice data" w p ls 1,\
     "section2d_simple.txt" u 1:2 i 1 t "interpolation" w l ls 2


set out
set terminal X
