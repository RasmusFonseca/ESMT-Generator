CXX = clang++

ESMTGenerator: main.cpp
	$(CXX) -o $@ $^

