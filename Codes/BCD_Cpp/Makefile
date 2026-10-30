CXX      = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra

bcd_rota: bcd_rota.cpp bcd_cetsp.cpp bcd_cetsp.hpp
	$(CXX) $(CXXFLAGS) bcd_cetsp.cpp bcd_rota.cpp -o bcd_rota

clean:
	rm -f bcd_rota points.txt

.PHONY: clean
