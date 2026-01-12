CXX=mpicxx
CPPFLAGS=-I/opt/local/include -I../include -I../external
LIBS=-lmpi
CXXFLAGS= -std=c++17 -Ofast -mtune=native -march=native
#CXXFLAGS= -std=c++17 -O2 -fsanitize=address -g
all : fput

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $<

fput : main.o $(patsubst %.cpp, %.o, $(wildcard *.cpp))
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $? -o $@ $(LIBS)

.PHONY : clean distclean

clean :
	$(RM) *.o

distclean : clean
	$(RM) fput
