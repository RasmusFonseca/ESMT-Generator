CXX = clang++

STPGenerator: main.cpp
	$(CXX) -o $@ $^

