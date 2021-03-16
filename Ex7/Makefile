args=10000 10 25 # serial / parallel / incl graham parallel / static scheduling average code time
#10000 10 25 -> 	~34s  /	5492ms	 / 5347ms/ 5378ms
#200 200 25 -> 		614ms /	478ms	 / 226ms / 202ms
#1000 10 25 -> 		370ms /	78ms  	 / 87ms  / 91ms
#1000 20 50 ->		1477ms/ 348ms	 / 294ms / 342ms

default: eigenvalues.exe
parallel: eigenvalues_parallel.exe
eigenvalues_parallel.exe: eigenvalues.cpp
	g++ -O3 -DNDEBUG -o eigenvalues_parallel.exe eigenvalues.cpp -fopenmp
eigenvalues.exe: eigenvalues.cpp
	g++ -O3 -DNDEBUG -o eigenvalues.exe eigenvalues.cpp
run: eigenvalues.exe eigenvalues_parallel.exe
	# ./eigenvalues.exe $(args)
	./eigenvalues_parallel.exe $(args)
all: default parallel run