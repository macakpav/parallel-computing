args=10000 10 25 # serial / parallel / incl vector product average code time
#10000 10 25 -> 	~34s  /	5492ms	 / 5492ms
#200 200 25 -> 		~600ms/	440ms	 / 
#1000 10 25 -> 		~367ms/	79ms  	 / 92ms
#1000 20 50 ->		1477ms/ 320ms	 / 310ms

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