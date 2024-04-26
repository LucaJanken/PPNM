# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRCS = main.cc ann.cc ../matrix/vector.cc ../matrix/matrix.cc ../minimization/min.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Default target
all: $(OUTPUT) ann_plot ann_derivative_plot ann_second_derivative_plot ann_integral_plot

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

# Rule for plotting network response against actual function
ann_plot: ann_results.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'ann_plot.svg'; \
	set title 'Network Response vs Actual Function'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot cos(5 * x - 1) * exp(-x * x) with lines linestyle 1 linecolor 'grey' title 'Actual Function', \
	     'ann_results.csv' using 1:2 with points pointtype 6 linecolor 'grey' title 'Training Data', \
	     'ann_results.csv' using 1:3 with lines linecolor 'red' title 'Network Response (n=5)'; \
"

# Rule for network response against derivative of actual function
ann_derivative_plot: ann_results.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'ann_derivative_plot.svg'; \
	set title 'Network Response vs Derivative of Actual Function'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot -2 * x * cos(5 * x - 1) * exp(-x * x) - 5 * sin(5 * x - 1) * exp(-x * x) with lines linestyle 1 linecolor 'grey' title 'Derivative of Actual Function', \
	     'ann_results.csv' using 1:4 with lines linecolor 'red' title 'Network Response (n=5)'; \
"
# Rule for network response against second derivative of actual function
ann_second_derivative_plot: ann_results.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'ann_second_derivative_plot.svg'; \
	set title 'Network Response vs Second Derivative of Actual Function'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot exp(-x * x) * ((4 * x * x - 27) * cos(5 * x - 1) + 20 * x * sin(5 * x - 1)) with lines linestyle 1 linecolor 'grey' title 'Second Derivative of Actual Function', \
	     'ann_results.csv' using 1:5 with lines linecolor 'red' title 'Network Response (n=5)'; \
"

# Extract the first value of column 8 from ann_results.csv
FIRST_VAL=$(shell awk 'NR==1 {print $$6}' ann_results.csv)

ann_integral_plot: ann_results.csv
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'ann_integral_plot.svg'; \
	set title 'Network Response vs Anti Derivative (Taylor Approximation)'; \
	set xlabel 'x'; \
	set ylabel 'y'; \
	set grid; \
	plot $(FIRST_VAL)-0.1877184794-0.045046858263+cos(1.0)*x+5.0/2.0*sin(1.0)*x**2-9.0/2.0*cos(1.0)*x**3-155.0/24.0*sin(1.0)*x**4+937.0/120.0*cos(1.0)*x**5+395.0/48.0*sin(1.0)*x**6-7799.0/1008.0*cos(1.0)*x**7-53215.0/8064.0*sin(1.0)*x**8+125087.0/24192.0*cos(1.0)*x**9+545549.0/145152.0*sin(1.0)*x**10-20393423.0/7983360.0*cos(1.0)*x**11-3478267.0/2128896.0*sin(1.0)*x**12+1231265381.0/1245404160.0*cos(1.0)*x**13+1982571053.0/3487131648.0*sin(1.0)*x**14-27192392077.0/87178291200.0*cos(1.0)*x**15-27417833143.0/167382319104.0*sin(1.0)*x**16+1174908885961.0/14227497123840.0*cos(1.0)*x**17+62190289289.0/1552090595328.0*sin(1.0)*x**18-91253890786099.0/4865804016353280.0*cos(1.0)*x**19-165135954461431.0/19463216065413120.0*sin(1.0)*x**20+844005190156393.0/227070854096486400.0*cos(1.0)*x**21 with lines linestyle 1 linecolor 'grey' title 'Anti Derivative (Taylor Approximation)', \
	     'ann_results.csv' using 1:6 with lines linecolor 'red' title 'Network Response (n=5)'; \
"

# Clean target
clean:
	rm -f $(TARGET) $(OBJS) $(OUTPUT) *.svg fit.log *.csv *.txt 

# Phony
.PHONY: all clean ann_plot
