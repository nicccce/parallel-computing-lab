CXX = g++
CXXFLAGS = -O3 -mavx2 -fopenmp -pthread
LIBS = -lz
TARGET = jaccard_cluster
SRC = wj.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) $(LIBS) -o $(TARGET)

clean:
	rm -f $(TARGET)
