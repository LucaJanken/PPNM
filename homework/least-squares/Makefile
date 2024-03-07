# Compiler settings
CXX = g++ 
CXXFLAGS = -std=c++11 -Wall -Wextra
LDLIBS = -lstdc++ -lm

# Source files
SRCS = main.cc ../matrix/matrix.cc ../matrix/vector.cc lsfit.cc ../linear_equations/QRGS.cc

# Datafile
DATAFILE = datafile.csv

# Object files
OBJS = $(SRC:.cc=.o)

# Target executable
TARGET = main

# Default targets
all: Out.txt datafit.svg

# Rule to generate the output file
Out.txt: $(TARGET)
	./$(TARGET) > $@

# Rule to create the target executable
$(TARGET): $(OBJS)
	$(CXX) $(SRCS) $(OBJS) $(LDLIBS) -o $(TARGET) $(CXXFLAGS)

# Rule to create the object files
%.o: %.cc
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Variables extracted from Out.txt
A = $(shell grep 'a =' Out.txt | awk '{print $$3}')
B = $(shell awk '/c:/{getline; print $$2}' Out.txt)
DA = $(shell grep 'Uncertainties:' Out.txt | awk '{print $$2}')
DB = $(shell grep 'Uncertainties:' Out.txt | awk '{print $$3}')


# Rule to generate the plot
datafit.svg: $(DATAFILE) Out.txt
	@echo "A = $(A), B = $(B), DA = $(DA), DB = $(DB)"
	gnuplot -e "\
	set terminal svg size 800,600 background 'white'; \
	set output 'datafit.svg'; \
	set title 'Decay of 224Ra with fit'; \
	set xlabel 't (days)'; \
	set ylabel 'Activity'; \
	f(x) = $(A)*exp(-$(B)*x); \
	fp(x) = ($(A)+$(DA))*exp(-($(B)-$(DB))*x); \
	fm(x) = ($(A)-$(DA))*exp(-($(B)+$(DB))*x); \
	plot '$(DATAFILE)' with errorbars title 'data', \
		 f(x) title sprintf('Fit: %.f*exp(-%.2f t)', $(A), $(B)), \
		 fp(x) with lines lt 1 lc 'red' title sprintf('Fit with plus sigma'), \
		 fm(x) with lines lt 1 lc 'green' title sprintf('Fit with minus sigma');"

# Clean
clean:
	rm -f $(TARGET) $(OBJS) Out.txt datafit.svg

# Phony
.PHONY: all clean


