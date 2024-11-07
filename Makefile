CXX=c++

CXXFLAGS= -std=c++17 -Ofast -mtune=native -march=native

all : fput

%.o : %.cpp
	$(CXX) -c $(CXXFLAGS) $<

fput : main.o $(patsubst %.cpp, %.o, $(wildcard *.cpp))
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $? -o $@ 

.PHONY : clean distclean

clean :
	$(RM) *.o

distclean : clean
	$(RM) flup