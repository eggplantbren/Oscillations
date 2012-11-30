# Development
#CFLAGS = -O0 -g -Wall -Wextra -ansi -pedantic

# Optimized
CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -DNDEBUG
LIBS = -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system

default:
	# Compile
	g++ $(CFLAGS) -c main.cpp
	g++ $(CFLAGS) -c FitSine.cpp
	g++ $(CFLAGS) -c Data.cpp

	# Link
	g++ -o main *.o $(LIBS)

clean:
	rm -f main *.o


