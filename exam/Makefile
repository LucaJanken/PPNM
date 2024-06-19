# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRCS = main.cc \
       RayleighQuotient.cc \
       EigenSolver.cc \
       HydrogenHamiltonian.cc \
       ScalingInvestigation.cc \
       ../homework/matrix/vector.cc \
       ../homework/matrix/matrix.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Default target
all: $(OUTPUT) plot

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

# Rule for plotting Scaling Results
plot: $(OUTPUT)
	@gnuplot -e "set terminal svg enhanced; \
	set output 'scaling.svg'; \
	set title 'Scaling Investigation'; \
	set xlabel 'Matrix Size'; \
	set ylabel 'Time (s)'; \
	set grid; \
	f(x) = a*x**2; \
	fit f(x) 'scaling_results.csv' using 1:2 via a; \
	plot 'scaling_results.csv' using 1:2 with linespoints title 'Evaluation time vs. Matrix Size', \
	f(x) title 'Fit: a*x^2';"

# Clean target
clean:
	rm -f $(TARGET) $(OBJS) $(OUTPUT) *.svg fit.log *.csv *.txt

# Phony
.PHONY: all clean plot
