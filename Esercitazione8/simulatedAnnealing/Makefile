CC = g++ -g
CFLAGS = -Wall -O0 --std=c++17

main.exe : main.o random.o params.o metropolis.o
	$(CC)  params.o random.o main.o metropolis.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
params.o : params.cpp
	$(CC) -c params.cpp -o params.o $(CFLAGS)
metropolis.o : metropolis.cpp
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
