Default: FORCE
	g++ -o prova main.cpp weisfeiler_lehman.cpp -lm -std=c++0x
debug: FORCE
	g++ -o prova main.cpp weisfeiler_lehman.cpp -lm -std=c++0x -g -DDEBUG -pthread
optimized: FORCE
	g++ -o prova main.cpp weisfeiler_lehman.cpp -lm -std=c++0x -O2
superoptimized: FORCE
	g++ -o prova main.cpp weisfeiler_lehman.cpp -lm -std=c++14 -O3 -pthread
debugoptimized: FORCE
	g++ -o prova main.cpp weisfeiler_lehman.cpp -lm -std=c++0x -g -DDEBUG -pthread -O3
FORCE:

