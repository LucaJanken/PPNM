# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Header file directory
INCLUDES = -I/home/lucajn/repos/PPNM/exercises/vec

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRC = main.cc /home/lucajn/repos/PPNM/exercises/vec/vec.cc

# Rule to generate the output file
$(OUTPUT): $(TARGET)						
	./$(TARGET) > $(OUTPUT)

# Rule to create the target executable
$(TARGET): $(SRC)
	$(CXX) $(SRC) $(INCLUDES) $(CXXFLAGS) $(LDLIBS) -o $(TARGET) 

# Clean target
clean:
	rm -f $(OUTPUT) $(TARGET)

# Phony targets
.PHONY: clean

