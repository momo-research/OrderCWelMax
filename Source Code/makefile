# Specify compiler and compilation options
CXX = g++
CXXFLAGS = -I./include -std=c++17 -Wall -O2 -fPIE -fopenmp
LDFLAGS = -fPIE -pie -fopenmp

# Executable file
TARGET = build/seqic

# All source files
SRCS = $(wildcard src/*.cpp)

# Generated object files (.o)
OBJS = $(patsubst src/%.cpp, build/%.o, $(SRCS))

# Default target: compile all files
all: $(TARGET)

# Link object files to generate the final executable
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@

# Compile source files to generate object files
build/%.o: src/%.cpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -rf build/*.o build/seqic