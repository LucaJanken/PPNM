# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRCS = main.cc roots.cc ../matrix/vector.cc ../matrix/matrix.cc ../linear_equations/QRGS.cc ../ode/rkode.cc ../splines/spline.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Default target
all: $(OUTPUT) wavefunction_plot convergence_plot convergence_improved_plot

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

# Plots
wavefunction_plot: wavefunction.csv
	gnuplot -e " \
	set terminal svg size 800,600 background 'white'; \
	set output 'wavefunction_plot.svg'; \
	set title 'Comparison of Numerical and Analytical Wavefunction for H (n = 0, s-wave)'; \
	set xlabel 'r (Bohr radius)'; \
	set ylabel 'f(r)'; \
	set grid; \
	plot 'wavefunction.csv' using 1:2 with points title 'Numerical Wavefunction', \
	x*exp(-x) with lines title 'Analytical Wavefunction';"

convergence_plot: convergence.csv
	gnuplot -e " \
	set terminal svg size 1200,800; \
	set output 'convergence_plot.svg'; \
	set multiplot layout 2,2 title 'Convergence Analysis'; \
	\
	set title 'rmin vs. Ermin'; \
	set xlabel 'rmin'; \
	set xrange [1:0]; \
	set ylabel 'Ermin'; \
	plot 'convergence.csv' using 1:2 with linespoints title 'Ermin', \
	-0.5 lt -1 lc 'red' lw 1 title 'Analytical'; \
	\
	set title 'rmax vs. Ermax'; \
	set xlabel 'rmax'; \
	set xrange [1:10]; \
	set ylabel 'Ermax'; \
	plot 'convergence.csv' using 3:4 with linespoints title 'Ermax', \
	-0.5 lt -1 lc 'red' lw 1 title 'Analytical'; \
	\
	set title 'acc vs. Eacc'; \
	set xlabel 'acc'; \
	set logscale x; \
	set xrange [500:1e-4]; \
	set ylabel 'Eacc'; \
	set yrange [-1:-0.4]; \
	plot 'convergence.csv' using 5:6 with linespoints title 'Eacc', \
	-0.5 lt -1 lc 'red' lw 1 title 'Analytical'; \
	\
	set title 'eps vs. Eeps'; \
	set xlabel 'eps'; \
	set xrange [500:1e-4]; \
	set ylabel 'Eeps'; \
	plot 'convergence.csv' using 7:8 with linespoints title 'Eeps', \
	-0.5 lt -1 lc 'red' lw 1 title 'Analytical'; \
	\
	unset multiplot;"

convergence_improved_plot: convergence_improved.csv
	gnuplot -e " \
	set terminal svg size 800,600; \
	set output 'convergence_improved_plot.svg'; \
	set title 'E vs rmax (Improved)'; \
	set xlabel 'rmax'; \
	set ylabel 'E'; \
	set yrange [-0.6:-0.4]; \
	set grid; \
	plot 'convergence_improved.csv' using 1:2 with linespoints title 'E (n = 1)';"
	

# Clean target
clean:
	rm -f $(TARGET) $(OBJS) $(OUTPUT) *.csv *.svg fit.log wavefunction_plot.svg convergence_plot.svg convergence_improved_plot.svg

# Phony
.PHONY: all clean wavefunction_plot