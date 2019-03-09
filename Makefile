output: main.o src/mesh.o
	g++ main.o mesh.o -o output
main.o:main.cc
	g++ -c main.cc -I ./Eigen -I ./include
src/mesh.o: src/mesh.cc
	g++ -c src/mesh.cc -I ./Eigen -I ./include

clean:
	rm *.o output
