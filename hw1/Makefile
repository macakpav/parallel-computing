default: 
	@echo "Only make gol, make clean and make all."
clean:
	@rm -f *.gol
all: clean gol_parallel run
	@python3 gol_visualization.py ????-??-??-??-??-??.gol
run: clean
	@ mpiexec -np 4 ./gol.exe 50 50 2 10
gol: gol_serial.cpp
	g++ -o gol.exe -O3 -std=c++1z gol_serial.cpp
gol_parallel: gol.cpp
	mpic++ -o gol.exe -g -Wall -std=c++17 gol.cpp
wipe: clean
	@rm -f *.exe