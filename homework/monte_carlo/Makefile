# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRCS = main.cc montecarlo.cc ../matrix/vector.cc ../matrix/matrix.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Default target
all: $(OUTPUT) circle_error_plot ellipse_error_plot circle_error_comparison_plot ellipse_error_comparison_plot

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
circle_error_plot: circle_data.csv
	gnuplot -e " \
	set terminal svg size 800,600 background 'white'; \
	set output 'circle_error_plot.svg'; \
	set title 'Circle area error'; \
	set xlabel 'Number of points'; \
	set xrange [100:25000]; \
	set ylabel 'Error'; \
	set yrange [0:0.5]; \
	set grid; \
	f(x) = a/sqrt(x); \
	fit f(x) 'circle_data.csv' using 1:4 via a; \
	plot 'circle_data.csv' using 1:3 with lines title 'Estimated error', \
	     'circle_data.csv' using 1:4 with lines title 'Actual error', \
	     f(x) with lines lw 2 title sprintf('Fitted 1/sqrt(N), a=%.2f', a);"

ellipse_error_plot: ellipse_data.csv
	gnuplot -e " \
	set terminal svg size 800,600 background 'white'; \
	set output 'ellipse_error_plot.svg'; \
	set title 'Ellipse area error'; \
	set xlabel 'Number of points'; \
	set xrange [100:25000]; \
	set ylabel 'Error'; \
	set yrange [0:0.5]; \
	set grid; \
	f(x) = a/sqrt(x); \
	fit f(x) 'ellipse_data.csv' using 1:4 via a; \
	plot 'ellipse_data.csv' using 1:3 with lines title 'Estimated error', \
	     'ellipse_data.csv' using 1:4 with lines title 'Actual error', \
	     f(x) with lines lw 2 title sprintf('Fitted 1/sqrt(N), a=%.2f', a);"

circle_error_comparison_plot: circle_data.csv circle_data_quasi.csv
	gnuplot -e " \
	set terminal svg size 800,600 background 'white'; \
	set output 'circle_error_comparison_plot.svg'; \
	set title 'Circle area error comparison'; \
	set xlabel 'Number of points (N)'; \
	set xrange [100:25000]; \
	set ylabel 'Error'; \
	set yrange [0:0.2]; \
	set grid; \
	plot 'circle_data.csv' using 1:4 with lines title 'Actual error (Pseudo-random)', \
	     'circle_data_quasi.csv' using 1:3 with lines title 'Estimated error (Quasi-random)', \
	     'circle_data_quasi.csv' using 1:4 with lines title 'Actual error (Quasi-random)';"

ellipse_error_comparison_plot: ellipse_data.csv ellipse_data_quasi.csv
	gnuplot -e " \
	set terminal svg size 800,600 background 'white'; \
	set output 'ellipse_error_comparison_plot.svg'; \
	set title 'Ellipse area error comparison'; \
	set xlabel 'Number of points (N)'; \
	set xrange [100:25000]; \
	set ylabel 'Error'; \
	set yrange [0:0.2]; \
	set grid; \
	plot 'ellipse_data.csv' using 1:4 with lines title 'Actual error (Pseudo-random)', \
	     'ellipse_data_quasi.csv' using 1:3 with lines title 'Estimated error (Quasi-random)', \
	     'ellipse_data_quasi.csv' using 1:4 with lines title 'Actual error (Quasi-random)';"

# Clean target
clean:
	rm -f $(OUTPUT) $(TARGET) $(OBJS) *.svg *.csv fit.log

# Phony targets
.PHONY: clean all circle_error_plot ellipse_error_plot circle_error_comparison_plot ellipse_error_comparison_plot