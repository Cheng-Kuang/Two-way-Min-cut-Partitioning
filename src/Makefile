CXX = g++
CXXFLAGS = -O3 -std=c++11 

SRC = partition.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = hw2
BIN_TARGET = ../bin/hw2

all: $(TARGET) $(BIN_TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ)

$(BIN_TARGET): $(OBJ)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $(BIN_TARGET) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(BIN_TARGET) *.o
