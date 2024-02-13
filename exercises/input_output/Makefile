# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Input file
INPUT = input.txt

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source file
SRC = main.cc

# Object file
OBJ = $(SRC:.cc=.o)

# Command-line argument numbers
CMD_LINE_NUMBERS = 1,2,3,4,5

# Standard input numbers
STDIN_NUMBERS = 5 6 7 8 9

all: $(TARGET) regenerate_output

# Combine tasks in a single output file
regenerate_output:
	# Task 1: Command-line arguments
	echo "# Task 1 (command-line arguments)" > $(OUTPUT)
	./$(TARGET) -numbers:$(CMD_LINE_NUMBERS) >> $(OUTPUT)
	# Task 2: Standard input
	echo "# Task 2 (standard input)" >> $(OUTPUT)
	echo -e "$(STDIN_NUMBERS)" | ./$(TARGET) >> $(OUTPUT)
	# Task 3: Input/output file
	echo "# Task 3 (input/output)" >> $(OUTPUT)
	./$(TARGET) -input:$(INPUT) -output:tmp_Out.txt
	cat tmp_Out.txt >> $(OUTPUT)
	rm -f tmp_Out.txt

# Rule to create the target executable
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) $(LDLIBS) -o $(TARGET)

# Rule for object file
$(OBJ): $(SRC)
	$(CXX) $(SRC) -c $(CXXFLAGS) -o $(OBJ)

# Clean target
clean:
	rm -f $(OUTPUT) $(TARGET) $(OBJ)

# Phony targets
.PHONY: clean all regenerate_output
