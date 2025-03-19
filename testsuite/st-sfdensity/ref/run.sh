cp predefines.h.1 predefines.h
make clean
make 1d

cp input.txt input.txt.run
echo "outprefix a0" >> input.txt.run
echo "temperature 0.0" >> input.txt.run
mpirun -np 20 ./st-wslda-1d input.txt.run

cp input.txt input.txt.run
echo "outprefix a2" >> input.txt.run
echo "temperature 0.2" >> input.txt.run
mpirun -np 20 ./st-wslda-1d input.txt.run

# ------------
cp predefines.h.2 predefines.h
make clean
make 1d

cp input.txt input.txt.run
echo "outprefix b0" >> input.txt.run
echo "temperature 0.0" >> input.txt.run
mpirun -np 20 ./st-wslda-1d input.txt.run

cp input.txt input.txt.run
echo "outprefix b2" >> input.txt.run
echo "temperature 0.2" >> input.txt.run
mpirun -np 20 ./st-wslda-1d input.txt.run

# ------------
cp predefines.h.3 predefines.h
make clean
make 1d

cp input.txt input.txt.run
echo "outprefix c0" >> input.txt.run
echo "temperature 0.0" >> input.txt.run
mpirun -np 20 ./st-wslda-1d input.txt.run

cp input.txt input.txt.run
echo "outprefix c2" >> input.txt.run
echo "temperature 0.2" >> input.txt.run
mpirun -np 20 ./st-wslda-1d input.txt.run
