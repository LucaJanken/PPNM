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
all: $(OUTPUT) plot plot_wavefunctions

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

# Plot targets
scaling.svg: scaling_results.csv $(OUTPUT)
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

wavefunctions.svg: wavefunctions.csv $(OUTPUT)
	@gnuplot -e "set terminal svg enhanced; \
	set output 'wavefunctions.svg'; \
	set title 'Hydrogen Wavefunctions'; \
	set xlabel 'r [arb. units]'; \
	set ylabel '|ψ(r)|^2'; \
	set xrange [0:12.3]; \
	set grid; \
	f1(x) = A * (2*(x-a)*exp(-(x-a)))**2; \
	f2(x) = B * (-(x-b)*(1-(x-b)/2)*exp(-(x-b)/2)/sqrt(2))**2; \
	fit f1(x) 'wavefunctions.csv' using 1:2 via A, a; \
	fit f2(x) 'wavefunctions.csv' using 1:3 via B, b; \
	plot 'wavefunctions.csv' using 1:2 with linespoints title 'Ground State', \
	     'wavefunctions.csv' using 1:3 with linespoints title 'First Excited State', \
	     f1(x) title 'Fitted Ground State (Analytical ψ_1)', \
	     f2(x) title 'Fitted First Excited State (Analytical ψ_2)';"

plot: scaling.svg

plot_wavefunctions: wavefunctions.svg

# Clean target
clean:
	rm -f $(TARGET) $(OBJS) $(OUTPUT) *.svg fit.log *.csv *.txt

# Phony
.PHONY: all clean plot plot_wavefunctions
