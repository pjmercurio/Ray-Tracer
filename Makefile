CC = g++
CFLAGS = -IEigen
CPPFLAGS = -I/opt/X11/include
LFLAGS = -L/opt/X11/lib
libs = -lX11
	
RM = /bin/rm -f 
all: main 

main: main.o 
	$(CC) $(CFLAGS) -o as2 main.o $(LFLAGS) ${libs}

example_01.o: example_01.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o
clean: 
	$(RM) *.o as2
 


