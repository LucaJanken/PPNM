# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output files
OUTPUT = Out.txt
PLOT_OUTPUT = benchmark.csv

# Target executable
TARGET = main

# Source files
SRCS = main.cc ../matrix/matrix.cc ../matrix/vector.cc QRGS.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Default target
all: $(OUTPUT) plot

# Rule to generate the output file and plot
$(OUTPUT): $(TARGET)
	./$(TARGET) > $(OUTPUT)

# Plot target
plot: $(PLOT_OUTPUT)
	gnuplot -e "\
	set terminal svg background "white"; \
	set terminal svg size 800,600; \
	set output 'benchmark_plot.svg'; \
	set title 'QR Decomposition Benchmark'; \
	set xlabel 'Matrix Size'; \
	set ylabel 'Time (seconds)'; \
	set grid; \
	f(x) = a*x**3; \
	fit f(x) '$(PLOT_OUTPUT)' using 1:2 via a; \
	plot '$(PLOT_OUTPUT)' using 1:2 with linespoints title 'Execution Time', \
	f(x) with lines title sprintf('Fit: %.3e x^3', a);\
	"

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
	rm -f $(OUTPUT) $(TARGET) $(OBJS) $(PLOT_OUTPUT) benchmark_plot.svg fit.log

# Phony targets
.PHONY: clean all plot
