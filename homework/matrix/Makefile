# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRCS = main.cc matrix.cc vector.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Rule to generate the output file
$(OUTPUT): $(TARGET)						
	./$(TARGET) > $(OUTPUT)

# Rule to create the target executable
$(TARGET): $(OBJS)	
	$(CXX) $(OBJS) $(LDLIBS) -o $(TARGET)

# Rule for main.o, depends on all header files
main.o: main.cc $(HDRS)
	$(CXX) $< -c $(CXXFLAGS) -o $@

# Generic rule for other object files, depends on their .cc and .h only
%.o: %.cc %.h
	$(CXX) $< -c $(CXXFLAGS) -o $@

# Clean target
clean:
	rm -f $(OUTPUT) $(TARGET) $(OBJS)

# Phony targets
.PHONY: clean