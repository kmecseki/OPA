# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra -pedantic -g

# Libraries
LIBS = -lgsl -lgslcblas -lfftw3 -static
INCLUDE_PATH = -I"C:/msys64/mingw64/include"
LIBRARY_PATH = -L"C:/msys64/mingw64/lib"

# Target executable
TARGET = OPAv6

# Source files
SRCS = run.cpp pulse.cpp utils.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default rule
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) $^ $(LIBRARY_PATH) $(LIBS) -o $@ 

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c $< -o $@

# Clean up generated files
clean:
	rm $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean


