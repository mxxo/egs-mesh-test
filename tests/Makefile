CXX      = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -g -O2 -I../

all: egs-mesh-tests

egs-mesh-tests: egs-mesh-tests.cpp ../msh_parser.h ../neighbour.h
		$(CXX) $(CXXFLAGS) egs-mesh-tests.cpp -o egs-mesh-tests
