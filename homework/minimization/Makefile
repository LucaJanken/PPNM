# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Output file
OUTPUT = Out.txt

# Target executable
TARGET = main

# Source files
SRCS = main.cc min.cc ../matrix/vector.cc ../matrix/matrix.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Header files
HDRS = $(patsubst %.cc,%.h,$(filter-out main.cc,$(SRCS)))

# Default target
all: $(OUTPUT) plot_higgs

# Rule to generate the output file by sending "higgs.data.csv" to the standard input
$(OUTPUT): $(TARGET) higgs.data.csv
	./$(TARGET) < higgs.data.csv > $(OUTPUT)

# Rule to create the target executable
$(TARGET): $(OBJS)	
	$(CXX) $(OBJS) $(LDLIBS) -o $(TARGET)

# Rule for main.o, depends on all header files
main.o: main.cc $(HDRS)
	$(CXX) $< -c $(CXXFLAGS) -o $@

# Generic rule for other object files, depends on their .cc and .h only
%.o: %.cc %.h
	$(CXX) $< -c $(CXXFLAGS) -o $@

# Rule for plotting Higgs data and fitting
plot_higgs: $(OUTPUT)
	@# Extract fitting parameters from line 23 and 26
	$(eval PARAMS := $(shell awk 'NR==23 {print $$1, $$2, $$3}' Out.txt))
	$(eval PARAMS2 := $(shell awk 'NR==26 {print $$1, $$2, $$3}' Out.txt))
	@# Use Gnuplot to plot
	@gnuplot -e "\
	set terminal svg enhanced; \
	set output 'higgs_fit.svg'; \
	set title 'Higgs Boson Signal with Breit-Wigner Fit'; \
	set xlabel 'Energy (GeV)'; \
	set ylabel 'Signal (arbitrary units)'; \
	set grid; \
	m = $(word 1, $(PARAMS)); \
	Gamma = $(word 2, $(PARAMS)); \
	A = $(word 3, $(PARAMS)); \
	m2 = $(word 1, $(PARAMS2)); \
	Gamma2 = $(word 2, $(PARAMS2)); \
	A2 = $(word 3, $(PARAMS2)); \
	f(x) = A / ((x - m)**2 + (Gamma / 2.0)**2); \
	f2(x) = A2 / ((x - m2)**2 + (Gamma2 / 2.0)**2); \
	set samples 1000; \
	plot 'higgs.data.csv' using 1:2:3 with errorbars title 'Experimental Data', \
	f(x) with lines title 'Breit-Wigner Fit (Symmetric Broyden)', \
	f2(x) with lines dashtype 2 linecolor 'red' title 'Breit-Wigner Fit (SR1)'; \
	"
# Clean target
clean:
	rm -f $(TARGET) $(OBJS) $(OUTPUT) *.svg fit.log

# Phony
.PHONY: all clean plot_higgs
