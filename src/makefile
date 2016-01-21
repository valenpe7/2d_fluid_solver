EXEC = fluid_solver
OBJS = main.o simulation_data.o cfd_simulation.o lodepng.o
SRCS = .
INCS = ./inc
LIBS = ./inc/lib

CC = mpicxx
CFLAGS = -Wall -O3 -c -std=c++0x
LFLAGS = -Wall -O3 -o $(EXEC)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

main.o: $(SRCS)/main.cpp $(INCS)/simulation_data.hpp $(INCS)/cfd_simulation.hpp $(LIBS)/lodepng/lodepng.h
	$(CC) $(CFLAGS) $(SRCS)/main.cpp

simulation_data.o: $(SRCS)/simulation_data.cpp $(INCS)/simulation_data.hpp
	$(CC) $(CFLAGS) $(SRCS)/simulation_data.cpp

cfd_simulation.o: $(SRCS)/cfd_simulation.cpp $(INCS)/cfd_simulation.hpp
	$(CC) $(CFLAGS) $(SRCS)/cfd_simulation.cpp

lodepng.o: $(LIBS)/lodepng/lodepng.cpp $(LIBS)/lodepng/lodepng.h
	$(CC) $(CFLAGS) $(LIBS)/lodepng/lodepng.cpp

clean:
	@rm *.o