# Development
#CFLAGS = -O0 -g -Wall -Wextra -ansi -pedantic

# Optimized
CFLAGS = -O2 -Wall -Wextra -ansi -pedantic -DNDEBUG
LIBS = -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system

default:
	# Compile
	g++ $(CFLAGS) -c *.cpp

	# Link
	g++ -o main *.o $(LIBS)
	rm -f main *.o

