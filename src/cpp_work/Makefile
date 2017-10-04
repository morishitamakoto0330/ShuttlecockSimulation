build:
	clang++ extraction.cpp -c -o extraction.o
	clang++ test.cpp -c -o test.o -std=c++11
	clang++ main.cpp -c -o main.o
	clang++ debug.cpp -c -o debug.o
	clang++ main.o extraction.o debug.o -o main `pkg-config --cflags opencv` `pkg-config --libs opencv`
	clang++ test.o extraction.o debug.o -o test `pkg-config --cflags opencv` `pkg-config --libs opencv`

run:
	./main

clean:
	#rm -f main test *.o *.out ../../res/image/test_*
	rm -f main test *.o *.out



