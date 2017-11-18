build:
	clang++ extraction.cpp -c -o extraction.o -std=c++11
	clang++ test.cpp -c -o test.o -std=c++11
	clang++ matrix.cpp -o matrix -std=c++11
	clang++ main.cpp -c -o main.o -std=c++11
	clang++ debug.cpp -c -o debug.o -std=c++11
	clang++ velocity_verlet.cpp -o velocity_verlet -std=c++11 `pkg-config --cflags --libs opencv`
	clang++ main.o extraction.o debug.o -o main `pkg-config --cflags opencv` `pkg-config --libs opencv`
	clang++ test.o extraction.o debug.o -o test `pkg-config --cflags opencv` `pkg-config --libs opencv`

run:
	./main

clean:
	#rm -f main test *.o *.out ../../res/image/test_*
	rm -f main test velocity_verlet matrix *.o *.out



