CC = g++ -g
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o salesman.o
	mpic++  salesman.o random.o main.o -o main.exe
main.o : main.cpp main.h
	mpic++ -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	
salesman.o : salesman.cpp salesman.h
	$(CC) -c salesman.cpp -o salesman.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
