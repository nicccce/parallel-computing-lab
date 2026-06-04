CXX = g++
CXXFLAGS = -O3 -mavx2 -fopenmp -pthread
EXTRA_CXXFLAGS ?= 
LIBS = -lz
TARGET ?= bin/jaccard_cluster
SRC = wj.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $(SRC) $(LIBS) -o $(TARGET)

clean:
	rm -f $(TARGET)
