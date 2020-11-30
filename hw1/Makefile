default: 
	@echo "Only make gol, make clean and make all."
clean:
	@rm -f *.gol
all: clean gol
	./gol.exe 50 50 2 10
	@python3 gol_visualization.py ????-??-??-??-??-??.gol
gol: gol_serial.cpp
	g++ -o gol.exe -O3 -std=c++1z gol_serial.cpp
gol_parallel: gol.cpp
	g++ -o gol.exe -g -Wall -Werror -std=c++1z gol.cpp
wipe: clean
	@rm -f *.exe