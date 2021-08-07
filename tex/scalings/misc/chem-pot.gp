# results for the uniform solution obtained on lattice 128x128x128  with DX=1

set term post enh eps color dashed 20
set output "chem-pot.eps"

set style line 1 pt 5 ps 2.0 lw 2. lt 1 lc 1   
set style line 2 pt 7 ps 2.0 lw 2. lt 1 lc 3
set style line 12 pt 2 ps 1.5 lw 4. lt 1 lc -1 dt 2 

set key top left font ",16"
set key spacing 1.5
set xrange[0:1]
set grid

set xlabel "k_F"
# set ylabel ""
set title "SLDA, uniform system, lattice 128^3 with DX=1"

a=0.01
b=0.397
linear(kF) = a*kF+b
fit linear(x) "chem-pot.txt" u 1:3 via a


plot "chem-pot.txt" u 1:2 t "E/E_{ffg}" w p ls 1, \
     "chem-pot.txt" u 1:3 t "{/Symbol m}/{/Symbol e}_{F}" w p ls 2,\
     0.397 not w l ls 1, linear(x) not w l ls 2

set out
set terminal X

system("convert -colorspace sRGB -density 200 chem-pot.eps chem-pot.png")
