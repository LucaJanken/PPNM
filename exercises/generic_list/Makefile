# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRCS = main.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = genlist.h

# Rule to generate the output file
$(OUTPUT): $(TARGET)
	./$(TARGET) < input.txt > $(OUTPUT)

# Rule to create the target executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDLIBS) -o $(TARGET)

# Rule for main.o, depends on all header files
main.o: main.cc $(HDRS)
	$(CXX) $< -c $(CXXFLAGS) -o $@

# Clean target
clean:
	rm -f $(OUTPUT) $(TARGET) $(OBJS)

# Phony targets
.PHONY: clean
