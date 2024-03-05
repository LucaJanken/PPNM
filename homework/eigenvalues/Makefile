# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -pthread
LDLIBS = -lstdc++ -lm

# Output files
OUTPUT = Out.txt
DELTA_R_OUTPUT = delta_r_vs_epsilon0.csv
RMAX_OUTPUT = rmax_vs_epsilon0.csv
EIGENFUNCTION_OUTPUT = eigenfunctions.csv
BENCHMARK_OUTPUT = benchmark.csv

# Target executable
TARGET = main
BENCHMARK_TARGET = benchmark

# Source files
SRCS = main.cc ../matrix/matrix.cc ../matrix/vector.cc JDWCS.cc
BENCHMARK_SRCS = benchmark.cc ../matrix/matrix.cc ../matrix/vector.cc JDWCS.cc

# Object files
OBJS = $(SRCS:.cc=.o)
BENCHMARK_OBJS = $(BENCHMARK_SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Default target
all: $(OUTPUT) delta_r_plot rmax_plot eigenfunction_plot run_benchmark benchmark_plot

# Rule to generate the output file and plot
$(OUTPUT): $(TARGET)
	./$(TARGET) > $(OUTPUT)

# Plot target
delta_r_plot: $(DELTA_R_OUTPUT)
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'delta_r_vs_epsilon0.svg'; \
	set title 'ε_0 vs Δr'; \
	set xlabel 'Δr'; \
	set ylabel 'ε_0'; \
	set grid; \
	plot '$(DELTA_R_OUTPUT)' using 1:2 with linespoints , \
	"

rmax_plot: $(RMAX_OUTPUT)
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'rmax_vs_epsilon0.svg'; \
	set title 'ε_0 vs r_{max}'; \
	set xrange [2:10]; \
	set xlabel 'r_{max}'; \
	set ylabel 'ε_0'; \
	set grid; \
	plot '$(RMAX_OUTPUT)' using 1:2 with linespoints , \
	"

eigenfunction_plot: $(EIGENFUNCTION_OUTPUT)
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'eigenfunctions.svg'; \
	set title 'Eigenfunctions'; \
	set xlabel 'r'; \
	set ylabel 'ψ'; \
	set grid; \
	set xrange [0:20]; \
	f1(x) = 2*x*exp(-x); \
	f2(x) = -x*(1-x/2)*exp(-x/2)/sqrt(2); \
	f3(x) = 2*x*(1-2*x/3 + 2*x*x/27)*exp(-x/3)/sqrt(27); \
	plot '$(EIGENFUNCTION_OUTPUT)' using 1:2 with points title 'ψ_1-data', \
	     '$(EIGENFUNCTION_OUTPUT)' using 1:3 with points title 'ψ_2-data', \
	     '$(EIGENFUNCTION_OUTPUT)' using 1:4 with points title 'ψ_3-data', \
	     f1(x) with lines title 'Analytical ψ_1' linewidth 3, \
	     f2(x) with lines title 'Analytical ψ_2' linewidth 3, \
	     f3(x) with lines title 'Analytical ψ_3' linewidth 3;"

benchmark_plot: $(BENCHMARK_OUTPUT)
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'benchmark.svg'; \
	set title 'JDWCS Benchmark (Parallelized)'; \
	set xlabel 'Matrix Size'; \
	set ylabel 'Time (seconds)'; \
	set grid; \
	f(x) = a*x**3; \
	fit f(x) '$(BENCHMARK_OUTPUT)' using 1:2 via a; \
	plot '$(BENCHMARK_OUTPUT)' using 1:2 with linespoints title 'Execution Time', \
	f(x) with lines title sprintf('Fit: %.3e x^3', a);\
	"

# Rule to create the target executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDLIBS) -o $(TARGET) $(CXXFLAGS)

# Rule for main.o, depends on all header files
main.o: main.cc $(HDRS)
	$(CXX) $< -c $(CXXFLAGS) -o $@

# Rule to generate the benchmark executable
$(BENCHMARK_TARGET): $(BENCHMARK_OBJS)
	$(CXX) $(BENCHMARK_OBJS) $(LDLIBS) -o $(BENCHMARK_TARGET) $(CXXFLAGS)

# Rule to run the benchmark
run_benchmark: $(BENCHMARK_TARGET)
	./$(BENCHMARK_TARGET)

# Generic rule for other object files, depends on their .cc and .h only
%.o: %.cc %.h
	$(CXX) $< -c $(CXXFLAGS) -o $@

# Clean target
clean:
	rm -f $(OUTPUT) $(TARGET) $(OBJS) $(BENCHMARK_TARGET) $(BENCHMARK_OBJS) $(DELTA_R_OUTPUT) $(RMAX_OUTPUT) $(EIGENFUNCTION_OUTPUT) $(BENCHMARK_OUTPUT) *.svg *.csv fit.log

# Phony targets
.PHONY: clean all delta_r_plot rmax_plot eigenfunction_plot run_benchmark benchmark_plot
