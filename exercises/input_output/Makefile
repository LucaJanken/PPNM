# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source file
SRC = main.cc

# Object file
OBJ = $(SRC:.cc=.o)

all: $(TARGET) regenerate_output

# Combine tasks in a single output file
regenerate_output:
	# Task 1: Command-line arguments
	echo "# Task 1 (command-line arguments)" > $(OUTPUT)
	./$(TARGET) -numbers:1,2,3,4,5 >> $(OUTPUT)
	# Task 2: Standard input
	echo "# Task 2 (standard input)" >> $(OUTPUT)
	echo -e "5\n6\n7\n8\n9" | ./$(TARGET) >> $(OUTPUT)
	# Task 3: Input/output file
	echo "# Task 3 (input/output)" >> $(OUTPUT)
	./$(TARGET) -input:input.txt -output:tmp_Out.txt
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
